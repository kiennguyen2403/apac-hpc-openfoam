#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:20:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N denorm_4n
#PBS -j oe

# Denormal diagnosis at 4 nodes, current Bonus setup, normal mapping, 15 steps per arm. Phase census
# (job 182140240) showed a per-subdomain efficiency factor (+-8% per cell/face in all loop phases) not explained
# by mesh counts or addressing locality. With -fp-model=precise FTZ/DAZ are off, so denormal operands take
# microcode assists. Two arms, both with tools/gamgCensus and per-rank perf counters (user mode):
#   base   - unchanged (tools/ftzDaz loaded with FTZ_DAZ_OFF=1, report only)
#   ftzdaz - DIAGNOSTIC ONLY: MXCSR FTZ+DAZ set at load (not IEEE; results differ slightly)
# perf: cycles, instructions, ASSISTS.FP, ASSISTS.ANY, MEM_LOAD_RETIRED.L3_MISS, FP_ARITH scalar double.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
RLIB=/scratch/ya41/kn3896/openfoam/tools/stepReuse/lib/linux64IcxDPInt32OptSPR/libstepReuse.so
SFF=/scratch/ya41/kn3896/openfoam/tools/simpleFoamFused/bin/linux64IcxDPInt32OptSPR/simpleFoamFused
CUR=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/diag/denormal_4node_${PBS_JOBID%%.*}
ZLIB=/scratch/ya41/kn3896/openfoam/tools/ftzDaz/libftzDaz.so
CLIB=/scratch/ya41/kn3896/openfoam/tools/gamgCensus/lib/linux64IcxDPInt32OptSPR/libgamgCensus.so
REFLOG=/scratch/ya41/kn3896/openfoam/contest-results/diag/vtune_ranks_4node_182139223/simpleFoam.log
STEPS=15
CASEDIR=$ROOT/case
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "DENORM_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CUR/logFiles/$l.log"; done
[[ -f $LIB && -f $FLIB && -f $RLIB && -f $CLIB && -f $ZLIB && -x $SFF && -f $REFLOG ]]
command -v perf
for l in "$LIB" "$FLIB" "$RLIB" "$CLIB" "$ZLIB" "$SFF"; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
cp "$CUR/system/fvSolution.fixedIter" fvSolution.basic
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.basic > fvSolution.bonus
[[ $(grep -c 'smoothSolverFixed;' fvSolution.bonus) == 3 && $(grep -c 'PCGFixed;' fvSolution.bonus) == 1 ]]
awk '{print} /^[[:space:]]*grad\(U\);/ && incache {print "    grad(p);"; print "    grad(k);"; print "    grad(omega);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.bonus > fvSolution.tmp && mv fvSolution.tmp fvSolution.bonus
[[ $(grep -c 'grad(p);' fvSolution.bonus) == 1 && $(grep -c 'grad(omega);' fvSolution.bonus) == 1 ]]
lscpu -p=CPU,CORE,SOCKET,NODE > lscpu.txt

mkdir -p "$CASEDIR/constant"
cp -a "$CUR/system" "$CASEDIR/"
cp fvSolution.bonus "$CASEDIR/system/fvSolution"
cp -a "$CUR"/constant/{transportProperties,turbulenceProperties} "$CASEDIR/constant/"
ln -s "$MESH" "$CASEDIR/constant/polyMesh"
for p in "$CUR"/processors416_*; do ln -s "$p" "$CASEDIR/${p##*/}"; done
cp -a "$CUR/dynamicCode" "$CASEDIR/dynamicCode"

set +eu
module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
    WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
set -eu
[[ $WM_OPTIONS == linux64IcxDPInt32OptSPR ]]
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FOAM_SIGFPE=true
export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx I_MPI_DEBUG=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_EXCLUDE_LIST I_MPI_PIN_DOMAIN_SIZE
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
PIN_NORMAL=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
PIN_PERM=$(tr ',' '\n' <<< "$PIN_NORMAL" | tac | paste -sd,)
[[ $(tr ',' '\n' <<< "$PIN_NORMAL" | wc -l) == 104 && $(tr ',' '\n' <<< "$PIN_PERM" | wc -l) == 104 ]]
HOSTS_NORMAL=$(sort -u "$PBS_NODEFILE" | paste -sd,)
HOSTS_PERM=$(sort -u "$PBS_NODEFILE" | tac | paste -sd,)
echo "hosts normal: $HOSTS_NORMAL"; echo "hosts perm:   $HOSTS_PERM"
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-4n-104ppn-scatter.dat
glibcT=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)

cat > "$ROOT/perf-wrap.sh" <<WRAP
#!/bin/bash
exec perf stat -x, -o "$ROOT/\$ARM/perf/rank\$PMI_RANK.txt" \
    -e cycles:u,instructions:u \
    -e cpu/event=0xc1,umask=0x02,name=assists_fp/u,cpu/event=0xc1,umask=0x1b,name=assists_any/u \
    -e cpu/event=0xd1,umask=0x20,name=l3_miss/u,cpu/event=0xc7,umask=0x01,name=fp_scalar_dp/u \
    -- env LD_PRELOAD="\$PRELOAD" "\$@"
WRAP
chmod +x "$ROOT/perf-wrap.sh"

cd "$CASEDIR"
cp "$TEAM/dicts/4node/controlDict" system/controlDict
sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
echo "libs (\"$LIB\");" >> system/controlDict

arm() {   # $1 name, $2 extra -genv args
    local name=$1; shift
    mkdir -p "$ROOT/$name/census" "$ROOT/$name/perf"
    cd "$CASEDIR"
    rm -rf postProcessing
    echo "ARM_START $name utc=$(date -u +%FT%TZ)"
    local status=0
    I_MPI_PIN_PROCESSOR_LIST=$PIN_NORMAL timeout --kill-after=30s 600s \
        mpirun -np 416 -ppn 104 -hosts "$HOSTS_NORMAL" "${glibcT[@]}" \
        -genv ARM "$name" -genv PRELOAD "$CLIB:$ZLIB:$FLIB:$RLIB" -genv FTZ_DAZ_REPORT 1 "$@" \
        -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1 -genv SIMPLEFOAM_FUSED_REPORT 1 \
        -genv STEP_REUSE_CENSUS "$ROOT/$name/census" \
        "$ROOT/perf-wrap.sh" "$SFF" -parallel "${fh[@]}" > "$ROOT/$name/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$ROOT/$name/"; fi
    [[ -z $(ls -d "$CASEDIR"/processors416_*/[1-9]* 2>/dev/null) ]]   # nothing written
    cd "$ROOT/$name"
    [[ $status == 0 ]]
    grep -qx End simpleFoam.log
    [[ $(grep -c '^Time = ' simpleFoam.log) == $STEPS ]]
    [[ $(ls census/rank*.txt | wc -l) == 416 && $(ls perf/rank*.txt | wc -l) == 416 ]]
    grep 'ftzDaz:' simpleFoam.log || true
    echo "RUN_SUCCESS arm=$name $(grep 'Average wall-clock time per time step' simpleFoam.log | tail -1)"
    cd "$ROOT"
}

arm base -genv FTZ_DAZ_OFF 1
pres() { awk -v n="$2" '/^Time = /{t++} t>n{exit} /Solving for p,/{print}' "$1"; }
cmp <(pres base/simpleFoam.log 10) <(pres "$REFLOG" 10) && echo "P_RESIDUALS_IDENTICAL arm=base steps=10"
for m in 'fusedVectorSolve: U fused path (v3 folded copy)' 'stepReuse: GaussSeidel: folded' \
    'stepReuse: laplacian: fused Gauss linear assembly' 'simpleFoamFused: pressure equation steps fused'; do
    grep -qF "$m" base/simpleFoam.log || { echo "MISSING_REPORT $m"; false; }
done
arm ftzdaz
paste <(pres base/simpleFoam.log $STEPS | sed 's/.*Initial residual = \([^,]*\),.*/\1/') \
      <(pres ftzdaz/simpleFoam.log $STEPS | sed 's/.*Initial residual = \([^,]*\),.*/\1/') > p_initial_residuals.txt
for a in base ftzdaz; do
    f=$(find $a/postProcessing -name 'coefficient.dat' | head -1)
    echo "CD_STEP$STEPS arm=$a $(grep '^# *Time' "$f" | tail -1 | tr -s ' \t' ' ' | cut -c1-80) | $(awk '!/^#/{l=$0} END{print l}' "$f" | cut -c1-120)"
done
echo "DENORM_COMPLETE root=$ROOT"
