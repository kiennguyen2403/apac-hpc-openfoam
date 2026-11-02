#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N gcc_skx_ab_4n
#PBS -j oe

# Compiler A/B at 4 nodes x 104 ranks, Basic tier (stock simpleFoam, unmodified v2512), all with Intel MPI 2021.17.2:
#   icx    ICX 2025.3.2  -O3 -xSAPPHIRERAPIDS -fp-model=precise  (our production build, platform linux64IcxDPInt32OptSPR)
#   gcc    GCC 14.2.0    -O3 -march=sapphirerapids              (build job, platform linux64GccDPInt32OptSPR)
#   gcc_skx GCC 14.2.0   -O3 -march=skylake-avx512            (GCC's equivalent of ICX -xCORE-AVX512, platform
#          linux64GccDPInt32OptSKX; needs GCCSKXOUT)
#   avx512 ICX 2025.3.2  -O3 -xCORE-AVX512 -fp-model=precise    (the bonus folder's build.pbs, platform
#          linux64IcxDPInt32Opt in apac-2026-basic/openfoam/src/OpenFOAM-v2512; arm skipped if that build is missing)
# Same case and runtime for every arm: scotch (4 104)/(1 0.5) partition (decomp-screen/case416_w05 through a linked copy),
# shipped fvSolution/fvSchemes, controlDict = shipped controlDict.noWrite with endTime STEPS, scatter pinning, tuned
# Allreduce file, glibc malloc settings, UCX setting, -fileHandler collated -ioRanks (0 104 208 312), -ppn 104
# (the Basic arm of scripts/final/final_runs.sh). One 3-step warmup per arm (loads/compiles the coded timer), then
# REPEATS (default 2) x arms of STEPS (default 20), fixed arm order as in ARMS (default "icx gcc gcc_skx").

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
: "${GCCOUT:?qsub -v GCCOUT=<build output dir with dynamicCode_gcc>}"
REPEATS=${REPEATS:-2}
STEPS=${STEPS:-20}
RES=/scratch/ya41/kn3896/openfoam/contest-results
REFD=$RES/compliance/reference
SHIP=$REFD/run/case/occDrivAerStaticMesh/system
K=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
BT=/scratch/ya41/apac-2026-basic/openfoam/src/OpenFOAM-v2512
CASE=$RES/decomp-screen/case416_w05
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
ROOT=$RES/compiler-ab/gcc_skx_ab_4node_${PBS_JOBID%%.*}
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "COMPILER_AB_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$REFD/run/4node/intel-scatter-cpus.awk" "$REFD/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for f in fvSolution.fixedIter fvSchemes controlDict.noWrite; do cmp "$SHIP/$f" "$CASE/system/$f"; done
cmp "$SHIP/fvSolution.fixedIter" "$CASE/system/fvSolution"

read -r -a ARMS <<< "${ARMS:-icx gcc gcc_skx}"
[[ -x $K/platforms/linux64GccDPInt32OptSPR/bin/simpleFoam ]]
[[ " ${ARMS[*]} " != *" gcc_skx "* ]] || { : "${GCCSKXOUT:?qsub -v GCCSKXOUT=<gcc skx build output>}"; [[ -x $K/platforms/linux64GccDPInt32OptSKX/bin/simpleFoam ]]; }
echo "ARMS ${ARMS[*]}"

for a in "${ARMS[@]}"; do
    d=$ROOT/case_$a
    mkdir -p "$d/constant"
    cp -a "$CASE/system" "$d/"
    sed -E "s/^(endTime\s+)[0-9]+;/\1$STEPS;/" "$SHIP/controlDict.noWrite" > "$d/system/controlDict"
    cp -a "$CASE"/constant/{transportProperties,turbulenceProperties} "$d/constant/"
    ln -s "$MESH" "$d/constant/polyMesh"
    for p in "$CASE"/processors416_*; do ln -s "$p" "$d/${p##*/}"; done
    case $a in
        icx) cp -a "$CASE/dynamicCode" "$d/" ;;
        gcc) cp -a "$GCCOUT/dynamicCode_gcc" "$d/dynamicCode" ;;
        gcc_skx) cp -a "$GCCSKXOUT/dynamicCode_gcc" "$d/dynamicCode" ;;
        avx512) ;;   # compiled by OpenFOAM in the warmup (this tree has its sources)
    esac
done

load_arm() {   # $1 arm; in a subshell
    set +eu
    module purge
    case $1 in
        icx|avx512) module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
                    export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx ;;
        gcc|gcc_skx) module load gcc/14.2.0 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
                    export CC=gcc CXX=g++ I_MPI_CC=gcc I_MPI_CXX=g++ ;;
    esac
    case $1 in
        icx)    source "$K/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR ;;
        gcc)    source "$K/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Gcc WM_CC=gcc WM_CXX=g++ WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR ;;
        gcc_skx) source "$K/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Gcc WM_CC=gcc WM_CXX=g++ WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSKX ;;
        avx512) source "$BT/etc/bashrc" ;;
    esac
    export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
    set -eu
    case $1 in
        icx) [[ $WM_OPTIONS == linux64IcxDPInt32OptSPR ]] ;;
        gcc) [[ $WM_OPTIONS == linux64GccDPInt32OptSPR ]] ;;
        gcc_skx) [[ $WM_OPTIONS == linux64GccDPInt32OptSKX ]] ;;
        avx512) [[ $WM_OPTIONS == linux64IcxDPInt32Opt ]] ;;
    esac
    export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FOAM_SIGFPE=true
    export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx I_MPI_DEBUG=0
    unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_EXCLUDE_LIST I_MPI_PIN_DOMAIN_SIZE
    export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
    export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
    I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
    export I_MPI_PIN_PROCESSOR_LIST
    unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
    export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-4n-104ppn-scatter.dat
}
glibcT=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)
ucxT=(-genv UCX_ZCOPY_THRESH inf -genv UCX_RNDV_THRESH inf)

run() {   # $1 arm, $2 label (warmup|r<n>)
    local a=$1 label=$2 dir=$ROOT/run_${1}_$2 steps=$STEPS
    mkdir "$dir"
    (
        load_arm "$a"
        cd "$ROOT/case_$a"
        if [[ $label == warmup ]]; then
            steps=3
            {
                echo "arm=$a wm_options=$WM_OPTIONS project=$(readlink -f "$WM_PROJECT_DIR")"
                echo "simpleFoam=$(command -v simpleFoam) sha256=$(sha256sum "$(command -v simpleFoam)" | cut -c1-16)"
                echo "libfiniteVolume sha256=$(sha256sum "$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib/libfiniteVolume.so" | cut -c1-16)"
                echo "compiler=$($CXX --version | head -1)"
                echo "cxxflags=$(wmake -show-cxxflags 2>/dev/null)"
                echo "mpi=$(mpirun --version 2>&1 | head -1)"
            } > "$dir/build-metadata.txt"
        fi
        sed -E "s/^(endTime\s+)[0-9]+;/\1$steps;/" "$SHIP/controlDict.noWrite" > system/controlDict
        echo "RUN_START arm=$a run=$label utc=$(date -u +%FT%TZ)"
        status=0
        timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 "${glibcT[@]}" "${ucxT[@]}" simpleFoam -parallel "${fh[@]}" \
            > "$dir/simpleFoam.log" 2>&1 || status=$?
        [[ -d postProcessing ]] && mv postProcessing "$dir/"
        [[ -z $(ls -d "$ROOT/case_$a"/processors416_*/[1-9]* 2>/dev/null) ]]
        [[ $status == 0 ]] || { echo "RUN_FAILED arm=$a run=$label status=$status"; tail -5 "$dir/simpleFoam.log"; exit 1; }
        grep -qx End "$dir/simpleFoam.log"
        [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $steps ]]
        grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log" || { echo "RUN_FAILED arm=$a run=$label: timer missing"; exit 1; }
        cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
        echo "RUN_SUCCESS arm=$a run=$label Cd_step$steps=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
    )
}

OK=()
for a in "${ARMS[@]}"; do if run "$a" warmup; then OK+=("$a"); else echo "ARM_DROPPED $a (warmup failed)"; fi; done
(( ${#OK[@]} >= 2 ))
for ((r = 1; r <= REPEATS; r++)); do for a in "${OK[@]}"; do run "$a" r$r || echo "RUN_FAILED_RECORDED $a r$r"; done; done

python3 - "$ROOT" "$REPEATS" "${OK[@]}" <<'PY'
import re, sys, glob, statistics as st
from pathlib import Path
root, reps, arms = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3:]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
def t(a, r):
    p = root / f'run_{a}_r{r}' / 'simpleFoam.log'
    m = pat.findall(p.read_text()) if p.exists() else []
    return float(m[-1]) if m else None
def coeffs(a):
    f = glob.glob(str(root / f'run_{a}_r1/postProcessing/forceCoeffsAll/*/coefficient*.dat'))
    return [l.split()[1] for l in open(f[0]) if l.strip() and not l.startswith('#')] if f else None
label = {'icx': 'ICX -xSAPPHIRERAPIDS', 'gcc': 'GCC -march=sapphirerapids', 'gcc_skx': 'GCC -march=skylake-avx512', 'avx512': 'ICX -xCORE-AVX512 (bonus folder)'}
T = {a: [x for x in (t(a, r) for r in range(1, reps + 1)) if x] for a in arms}
b = st.mean(T['icx'])
C = {a: coeffs(a) for a in arms}
out = ['# Compiler A/B, Basic tier (stock simpleFoam), Intel MPI 2021.17.2, 4 nodes x 104, scotch, 20 steps', '',
       '| Build | repeats (s/step) | mean | vs ICX-SPR | Cd step 20 | coefficient.dat vs ICX-SPR |', '|---|---|---:|---:|---|---|']
for a in arms:
    if not T[a]:
        out.append(f'| {label[a]} | failed | | | | |'); continue
    same = 'identical' if C[a] == C['icx'] else (f'differs (max {max(abs(float(x)/float(y)-1) for x, y in zip(C[a], C["icx"]))*100:.4f}%)' if C[a] and C['icx'] else '?')
    out.append(f'| {label[a]} | {" / ".join(f"{x:.4f}" for x in T[a])} | {st.mean(T[a]):.4f} | {100*(st.mean(T[a])/b-1):+.2f}% | {C[a][-1] if C[a] else "?"} | {same if a != "icx" else "-"} |')
(root / 'REPORT.md').write_text('\n'.join(out) + '\n'); print('\n'.join(out))
PY
python3 /scratch/ya41/kn3896/openfoam/tools/perfModel/ab_steps.py "$ROOT" 2>/dev/null | tail -8 || true
for a in "${OK[@]}"; do d=case_$a; n=$(find "$d" | wc -l); tar czf "$d.tar.gz" "$d" && [[ $(tar tzf "$d.tar.gz" | wc -l) == "$n" ]] && rm -rf "$d"; done
echo "COMPILER_AB_COMPLETE root=$ROOT"
