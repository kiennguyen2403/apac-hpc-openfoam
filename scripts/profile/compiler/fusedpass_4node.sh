#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N fusedpass_4n
#PBS -j oe

# Source change for the memory-bound solver (BONUS tier: code change): with fixed-iteration
# controls (tolerance 0, relTol 0), smoothSolver computes a full residual (matrix pass + field
# allocation + global reduction) after every sweep only to test a convergence that cannot happen;
# PCG computes a residual norm (global reduction) every iteration. tools/fixedIterSolvers adds
# smoothSolverFixed / PCGFixed which compute these once after the last sweep/iteration. Same
# iterates, so results must be bit-identical.
# Paired test, 4 nodes x 104 ranks, prepared 4-node case416_w05 (scotch (4 104)/(1 0.5), no prep), 10 steps per arm, current best
# setup (ICX OptSPR + Intel MPI, scatter pinning, tuned Allreduce, glibcT).
#   base  : team fvSolution (smoothSolver for U/k/omega, PCG on the GAMG coarsest level)
#   fixed : same with smoothSolverFixed / PCGFixed, library loaded via controlDict libs

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
CASE=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
FLIB1=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib_v1/libfusedVectorSolve.so
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/fusedu/passes_4node_${PBS_JOBID%%.*}
STEPS=20
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "FUSEDPASS_FAILED line=$LINENO status=$?"; cp "$ROOT/fvSolution.orig" "$CASE/system/fvSolution" 2>/dev/null || true' ERR
cd "$ROOT"
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/logFiles/$l.log"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
ls "$CASE/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
[[ -f $LIB && -f $FLIB && -f $FLIB1 ]]
cp "$CASE/system/fvSolution" fvSolution.orig
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.orig > fvSolution.fixed
echo "fvSolution changes (base -> fixed):"; diff fvSolution.orig fvSolution.fixed || true
[[ $(grep -c 'smoothSolverFixed;' fvSolution.fixed) == 3 && $(grep -c 'PCGFixed;' fvSolution.fixed) == 1 ]]
if [[ -d $CASE/postProcessing ]]; then mv "$CASE/postProcessing" "$ROOT/postProcessing.before"; fi

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
I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
export I_MPI_PIN_PROCESSOR_LIST
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-4n-104ppn-scatter.dat
glibcT=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)

[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
mpirun -np 416 -ppn 104 hostname | sort | uniq -c > placement.txt
if [[ $(wc -l < placement.txt) != 4 ]] || awk '$1 != 104 {bad=1} END {exit !bad}' placement.txt; then
    echo "PLACEMENT_WRONG:"; cat placement.txt; exit 1
fi

run() {   # $1 arm (base|fixed)
    local a=$1 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    cd "$CASE"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    cp "$ROOT/fvSolution.fixed" system/fvSolution
    echo "libs (\"$LIB\");" >> system/controlDict
    local pre=()
    if [[ $a == v1 ]]; then pre=(-genv LD_PRELOAD "$FLIB1" -genv FUSED_VECTOR_SOLVE_REPORT 1); else pre=(-genv LD_PRELOAD "$FLIB" -genv FUSED_VECTOR_SOLVE_REPORT 1); fi
    echo "RUN_START arm=$a utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 "${glibcT[@]}" "${pre[@]}" simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cp "$ROOT/fvSolution.orig" system/fvSolution
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "ACTIVE arm=$a: $(grep -m1 fusedVectorSolve "$dir/simpleFoam.log" || echo none)"
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
run v1
run v2
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"   # case restored

if diff <(grep -E 'Solving for|continuity errors' run_v1/simpleFoam.log) <(grep -E 'Solving for|continuity errors' run_v2/simpleFoam.log) > solver.diff; then
    echo "CHECK_SOLVER_LINES_IDENTICAL: yes ($(grep -c 'Solving for' run_v1/simpleFoam.log) lines)"
else
    echo "CHECK_SOLVER_LINES_IDENTICAL: NO"; head -6 solver.diff
fi
if cmp -s <(grep -v '^#' run_v1/postProcessing/forceCoeffsAll/*/coefficient*.dat) <(grep -v '^#' run_v2/postProcessing/forceCoeffsAll/*/coefficient*.dat); then
    echo "CHECK_FORCE_COEFFS_IDENTICAL: yes"
else
    echo "CHECK_FORCE_COEFFS_IDENTICAL: NO"
fi
python3 - "$ROOT" "$STEPS" <<'PY'
import re, sys
from pathlib import Path
root, steps = Path(sys.argv[1]), sys.argv[2]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
t = {a: float(pat.findall((root / f'run_{a}' / 'simpleFoam.log').read_text())[-1]) for a in ('v1', 'v2')}
lines = [f'# fusedVectorSolve v2 (initial/final residual passes fused) vs v1, plugin on, 4 nodes x 104, occDrivAerStaticMesh 65M, {steps} steps, same job', '',
         '| Arm | s/step | vs v1 |', '|---|---:|---:|',
         f"| v1 (fused sweeps) | {t['v1']:.4f} | - |",
         f"| v2 (+ fused residual/norm passes) | {t['v2']:.4f} | {100 * (t['v2'] / t['v1'] - 1):+.2f}% |"]
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n'); print('\n'.join(lines))
PY
echo "FUSEDPASS_COMPLETE root=$ROOT"
