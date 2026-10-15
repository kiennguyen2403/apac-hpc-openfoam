#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N bench_motorbike_spr
#PBS -j oe

# One-node compiler comparison on HPC_motorbike Medium (17.4M cells, 104-way decomposition).
# Builds: vanilla (module OpenFOAM/2512: GCC 8.5, OpenMPI 4.1.7) and the custom
# v2512 OptSPR builds gcc / clang / icx (Intel MPI 2021.17.2).
# 40 simpleFoam iterations per run; metric = contest timer
# "Average wall-clock time per time step" (from the second step on).
# Per build: warmup (2 steps, compiles the coded timer), two timed repeats in
# mirrored order, then one perf trial (perf stat on all ranks except the
# perf record ranks, perf record on ranks 0 26 52 78). Timed runs are perf-free.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
[[ ${PBS_NCPUS:-0} == 104 ]]

MESH=/scratch/ya41/kn3896/openfoam/contest-results/compiler-comparison/HPC_motorbike
TIMER=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/system/FOwallClockTimeStatistics
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
SCRIPTS=/scratch/ya41/kn3896/openfoam/scripts/profile/compiler
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/compiler-comparison/bench_motorbike_${PBS_JOBID%%.*}
ITERS=40
NP=104
VARIANTS=(vanilla gcc clang icx)
RECORD_RANKS="0 26 52 78"

grep -q '^MESH_SUCCESS cells=17444758 processors=104' /scratch/ya41/kn3896/openfoam/logs/mesh/mesh_motorbike_M.o180709841
mkdir "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "BENCHMARK_FAILED line=$LINENO status=$?"' ERR
sha256sum "$SCRIPTS/bench_motorbike_spr.sh" "$SCRIPTS/bench_motorbike_summary.py" > "$ROOT/workflow.sha256"
lscpu > "$ROOT/lscpu.txt"
sort -u "$PBS_NODEFILE" > "$ROOT/hosts"
cat /proc/sys/kernel/perf_event_paranoid > "$ROOT/perf_event_paranoid.txt"

# --- environments -----------------------------------------------------------
# Called in a subshell; leaves simpleFoam etc. and an mpirun on PATH and sets MPI_RUN.
load_variant() {
    set +eu
    module purge
    unset FOAM_EXTRA_CFLAGS FOAM_EXTRA_CXXFLAGS FOAM_BUILDROOT LD_PRELOAD
    export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
    local comp cc cxx flags mods
    case $1 in
        vanilla)
            module load OpenFOAM/2512
            # Site bashrc is pre-configured; the coded timer is compiled with system GCC.
            export WM_COMPILER_TYPE=system WM_CC=gcc WM_CXX=g++ CC=gcc CXX=g++
            MPI_RUN=(mpirun -np $NP --map-by core --bind-to core)
            set -eu
            return 0
            ;;
        gcc)   mods="gcc/14.2.0 binutils/2.44";      comp=Gcc;   cc=gcc;   cxx=g++;     flags="-march=sapphirerapids" ;;
        clang) mods="clang/19.1.6";                  comp=Clang; cc=clang; cxx=clang++; flags="-march=sapphirerapids -ffp-model=precise" ;;
        icx)   mods="intel-compiler-llvm/2025.3.2";  comp=Icx;   cc=icx;   cxx=icpx;    flags="-xSAPPHIRERAPIDS -fp-model=precise" ;;
    esac
    module load $mods intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
    export CC=$cc CXX=$cxx I_MPI_CC=$cc I_MPI_CXX=$cxx
    export FOAM_EXTRA_CFLAGS="$flags" FOAM_EXTRA_CXXFLAGS="$flags"
    source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=$comp WM_CC=$cc WM_CXX=$cxx \
        WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
    export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
    export FOAM_EXTRA_CFLAGS="$flags" FOAM_EXTRA_CXXFLAGS="$flags"
    export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm
    # Same placement as OpenMPI --map-by core --bind-to core: rank i on core i.
    export I_MPI_PIN=1 I_MPI_PIN_PROCESSOR_LIST=0-$((NP-1))
    unset I_MPI_PIN_DOMAIN I_MPI_PIN_ORDER
    MPI_RUN=(mpirun -np $NP)
    set -eu
    [[ $WM_OPTIONS == linux64${comp}DPInt32OptSPR ]]
}

record_env() {   # $1 = output dir
    {
        echo "variant=$v wm_options=$WM_OPTIONS"
        echo "simpleFoam=$(command -v simpleFoam)"
        echo "mpirun=$(command -v mpirun)"
        echo "mpi=$(mpirun --version 2>&1 | head -1)"
        echo "compiler=$($CXX --version 2>&1 | head -1)"
        echo "foam_extra_cxxflags=${FOAM_EXTRA_CXXFLAGS:-}"
        readelf -p .comment "$(command -v simpleFoam)" 2>/dev/null | grep -v '^$'
    } > "$1/build-metadata.txt"
    ldd "$(command -v simpleFoam)" > "$1/simpleFoam-ldd.txt"
    ! grep -q 'not found' "$1/simpleFoam-ldd.txt"
}

# --- shared prepared case (vanilla build, once) -----------------------------
PREP=$ROOT/prep
mkdir "$PREP"
cp -a "$MESH/system" "$PREP/"
mkdir "$PREP/constant"
for f in "$MESH"/constant/*; do
    if [[ ${f##*/} == polyMesh ]]; then ln -s "$f" "$PREP/constant/polyMesh"; else cp -a "$f" "$PREP/constant/"; fi
done
for p in "$MESH"/processor*; do cp -a "$p" "$PREP/"; done
cp "$TIMER" "$PREP/system/FOwallClockTimeStatistics"
(
    v=vanilla; load_variant vanilla
    cd "$PREP"
    for kv in "startFrom startTime" "startTime 0" "endTime $ITERS" "writeInterval 100000" "runTimeModifiable false"; do
        foamDictionary -entry ${kv% *} -set "${kv#* }" system/controlDict > /dev/null
    done
    ! grep -q '^functions' system/controlDict
    printf '\nfunctions\n{\n    #include "FOwallClockTimeStatistics"\n}\n' >> system/controlDict
    foamDictionary -entry functions -keywords system/controlDict | grep -qx wallClockTimeStatistics
    "${MPI_RUN[@]}" renumberMesh -constant -overwrite -parallel > log.renumberMesh 2>&1
    "${MPI_RUN[@]}" potentialFoam -writephi -parallel > log.potentialFoam 2>&1
    grep -q '^End' log.renumberMesh && grep -q '^End' log.potentialFoam
)
[[ $(find "$PREP" -maxdepth 1 -type d -name 'processor[0-9]*' | wc -l) == $NP ]]
(cd "$PREP" && find processor*/0 -type f | sort | xargs sha256sum) > "$PREP/initial-fields.sha256"
echo "PREPARED renumberMesh+potentialFoam (vanilla)"

# Per-variant case: own system/ (dynamicCode), read-only shared mesh, copied initial fields.
for v in "${VARIANTS[@]}"; do
    c=$ROOT/case_$v
    mkdir -p "$c"
    cp -a "$PREP/system" "$c/"
    ln -s "$PREP/constant" "$c/constant"
    for ((r=0; r<NP; r++)); do
        mkdir "$c/processor$r"
        ln -s "$PREP/processor$r/constant" "$c/processor$r/constant"
        cp -a "$PREP/processor$r/0" "$c/processor$r/"
    done
done

check_fields() {   # fields must be untouched by noWrite runs
    (cd "$ROOT/case_$1" && find processor*/0 -type f | sort | xargs sha256sum) | cmp - "$PREP/initial-fields.sha256"
}

# --- runs ----------------------------------------------------------------------
cat > "$ROOT/perf-wrap.sh" <<'EOF'
#!/bin/bash
r=${OMPI_COMM_WORLD_RANK:-${PMI_RANK:?no MPI rank variable}}
case " $RECORD_RANKS " in
    *" $r "*) exec perf record -q -F 499 -e cycles:u -o "$PERF_DIR/record.$r.data" -- "$@" ;;
    *) exec perf stat -x, -o "$PERF_DIR/stat.$r.csv" \
           -e cycles:u,instructions:u,branch-misses:u,cache-misses:u \
           -e cpu/event=0xc7,umask=0x01,name=fp_scalar_double/u \
           -e cpu/event=0xc7,umask=0x04,name=fp_128b_packed_double/u \
           -e cpu/event=0xc7,umask=0x10,name=fp_256b_packed_double/u \
           -e cpu/event=0xc7,umask=0x40,name=fp_512b_packed_double/u \
           -- "$@" ;;
esac
EOF
chmod +x "$ROOT/perf-wrap.sh"

run() {   # $1 variant, $2 label (warmup|t1|t2|perf)
    local v=$1 label=$2 dir=$ROOT/run_${1}_$2
    mkdir "$dir"
    (
        load_variant "$v"
        cd "$ROOT/case_$v"
        check_fields "$v"
        local steps=$ITERS wrap=()
        if [[ $label == warmup ]]; then
            steps=2
            record_env "$dir"
            export I_MPI_DEBUG=4   # rank pinning report (Intel MPI)
            [[ $v == vanilla ]] && MPI_RUN+=(--report-bindings)
            # perf must work inside the job before the perf trials are trusted.
            perf stat -x, -e instructions:u -o "$dir/perf-selftest.csv" true
        fi
        foamDictionary -entry endTime -set $steps system/controlDict > /dev/null
        if [[ $label == perf ]]; then
            export PERF_DIR=$dir RECORD_RANKS
            MPI_RUN+=("$ROOT/perf-wrap.sh")
        fi
        echo "RUN_START variant=$v run=$label utc=$(date -u +%FT%TZ)"
        local t0=$SECONDS
        timeout 1500 "${MPI_RUN[@]}" simpleFoam -parallel > "$dir/simpleFoam.log" 2>&1
        echo "elapsed_s=$((SECONDS-t0))" > "$dir/elapsed.txt"
        foamDictionary -entry endTime -set $ITERS system/controlDict > /dev/null
        grep -qx End "$dir/simpleFoam.log"
        [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $steps ]]
        [[ -d postProcessing ]] && mv postProcessing "$dir/"
        check_fields "$v"
        if [[ $label == perf ]]; then
            for r in $RECORD_RANKS; do
                perf report -i "$dir/record.$r.data" --stdio --sort dso,sym --percent-limit 0.3 \
                    > "$dir/report.$r.txt" 2>/dev/null
            done
        fi
        echo "RUN_SUCCESS variant=$v run=$label $(grep -m1 'Average wall-clock time per time step' "$dir/simpleFoam.log" || true)"
    )
}

for v in "${VARIANTS[@]}"; do run "$v" warmup; done
for v in vanilla gcc clang icx; do run "$v" t1; done
for v in icx clang gcc vanilla; do run "$v" t2; done
for v in "${VARIANTS[@]}"; do run "$v" perf; done

module purge
python3 "$SCRIPTS/bench_motorbike_summary.py" "$ROOT" "$RECORD_RANKS"
sha256sum -c --strict "$ROOT/workflow.sha256" > /dev/null
echo "BENCHMARK_COMPLETE root=$ROOT"
