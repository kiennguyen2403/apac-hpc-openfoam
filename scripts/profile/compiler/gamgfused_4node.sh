#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N gamgfused_4n
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
GLIB=/scratch/ya41/kn3896/openfoam/tools/gamgFused/lib/linux64IcxDPInt32OptSPR/libgamgFused.so
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/gamgfused/bench_4node_${PBS_JOBID%%.*}
STEPS=20
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "GAMGFUSED_FAILED line=$LINENO status=$?"; cp "$ROOT/fvSolution.orig" "$CASE/system/fvSolution" 2>/dev/null || true' ERR
cd "$ROOT"
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/logFiles/$l.log"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
ls "$CASE/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
[[ -f $LIB && -f $FLIB && -f $GLIB ]]
cp "$CASE/system/fvSolution" fvSolution.orig
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.orig > fvSolution.fixed
echo "fvSolution changes (base -> fixed):"; diff fvSolution.orig fvSolution.fixed || true
[[ $(grep -c 'smoothSolverFixed;' fvSolution.fixed) == 3 && $(grep -c 'PCGFixed;' fvSolution.fixed) == 1 ]]
# grad(p) added to the existing cache block (which already holds grad(U))
awk '{print} /^[[:space:]]*grad\(U\);/ && incache {print "    grad(p);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.fixed > fvSolution.gradp
[[ $(grep -c 'grad(p);' fvSolution.gradp) == 1 ]]
echo "fvSolution change (fixed -> gradp):"; diff fvSolution.fixed fvSolution.gradp || true
sed -e 's/^\(\s*solver\s\+\)GAMG;/\1GAMGFused;/' fvSolution.gradp > fvSolution.gamg
[[ $(grep -c 'GAMGFused;' fvSolution.gamg) == 1 ]]
echo "fvSolution change (gradp -> gamg):"; diff fvSolution.gradp fvSolution.gamg || true
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
    if [[ $a == gamg ]]; then cp "$ROOT/fvSolution.gamg" system/fvSolution; else cp "$ROOT/fvSolution.gradp" system/fvSolution; fi
    if [[ $a == gamg ]]; then echo "libs (\"$LIB\" \"$GLIB\");" >> system/controlDict; else echo "libs (\"$LIB\");" >> system/controlDict; fi
    local pre=()
    pre=(-genv LD_PRELOAD "$FLIB" -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv GAMG_FUSED_REPORT 1)
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
    echo "ACTIVE arm=$a: $(grep -E 'fusedVectorSolve|gamgFused' "$dir/simpleFoam.log" | sort -u | tr '\n' ' ')"
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
run base
run gamg
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"   # case restored

if diff <(grep -E 'Solving for|continuity errors' run_base/simpleFoam.log) <(grep -E 'Solving for|continuity errors' run_gamg/simpleFoam.log) > solver.diff; then
    echo "CHECK_SOLVER_LINES_IDENTICAL: yes ($(grep -c 'Solving for' run_base/simpleFoam.log) lines)"
else
    echo "CHECK_SOLVER_LINES_IDENTICAL: NO"; head -6 solver.diff
fi
if cmp -s <(grep -v '^#' run_base/postProcessing/forceCoeffsAll/*/coefficient*.dat) <(grep -v '^#' run_gamg/postProcessing/forceCoeffsAll/*/coefficient*.dat); then
    echo "CHECK_FORCE_COEFFS_IDENTICAL: yes"
else
    echo "CHECK_FORCE_COEFFS_IDENTICAL: NO"
fi
python3 - "$ROOT" "$STEPS" <<'PY'
import re, sys
from pathlib import Path
root, steps = Path(sys.argv[1]), sys.argv[2]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
t = {a: float(pat.findall((root / f'run_{a}' / 'simpleFoam.log').read_text())[-1]) for a in ('base', 'gamg')}
lines = [f'# GAMGFused (fused residual/restriction/scaling, fixed-iteration norm) vs GAMG, Bonus setup incl. grad(p) cache, 4 nodes x 104, occDrivAerStaticMesh 65M, {steps} steps, same job', '',
         '| Arm | s/step | vs base |', '|---|---:|---:|',
         f"| base (plugin + fused v2 + grad(p) cache) | {t['base']:.4f} | - |",
         f"| gamg (+ GAMGFused for p) | {t['gamg']:.4f} | {100 * (t['gamg'] / t['base'] - 1):+.2f}% |"]
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n'); print('\n'.join(lines))
PY
echo "GAMGFUSED_COMPLETE root=$ROOT"
