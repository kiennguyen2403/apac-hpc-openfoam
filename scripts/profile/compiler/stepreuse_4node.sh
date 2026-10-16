#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:35:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N stepreuse_4n
#PBS -j oe

# stepReuse bundle (BONUS tier: code change, all bit-identical):
#   - GaussSeidel: per-sweep bPrime = source copy folded into the previous sweep (p via the
#     preloaded smoother, U in fusedVectorSolve v3, k/omega via smoothSolverFixed v2 which does
#     all sweeps in one smoother call)
#   - kOmegaSST correct(): flux divergence once instead of 3x, upwind convection coefficients of
#     omega/k built once
#   - fvSolution cache { grad(k); grad(omega); } (stock code): the gradients computed for CDkOmega
#     are reused by the corrected laplacian terms instead of being computed again
#   - wallShearStress: evaluated in write() only (never written in the scored run)
# Paired test, 4 nodes x 104 ranks, prepared case416_w05, 20 steps per arm, same nodes:
#   base     : current Bonus (fixedIterSolvers v1, fusedVectorSolve v2, grad(p) cache)
#   cachekw  : base + cache grad(k), grad(omega) (no code change)
#   all      : fixedIterSolvers v2, fusedVectorSolve v3, LD_PRELOAD libstepReuse, + grad(k), grad(omega) cache
#   nokomega : all with STEP_REUSE_OFF=komega (attribution, only with -v ARMS="... nokomega")
#   nogs     : all with STEP_REUSE_OFF=gs (attribution: no folded copy in GaussSeidel or fused U, only with -v ARMS)

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
CASE=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
T=/scratch/ya41/kn3896/openfoam/tools
LIB1=$T/fixedIterSolvers/lib_v1/libfixedIterSolvers.so
LIB2=$T/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB2=$T/fusedVectorSolve/lib_v2/libfusedVectorSolve.so
FLIB3=$T/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
RLIB=$T/stepReuse/lib/linux64IcxDPInt32OptSPR/libstepReuse.so
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/stepreuse/bench_4node_${PBS_JOBID%%.*}
STEPS=20
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "STEPREUSE_FAILED line=$LINENO status=$?"; cp "$ROOT/fvSolution.orig" "$CASE/system/fvSolution" 2>/dev/null || true' ERR
cd "$ROOT"
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/logFiles/$l.log"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
ls "$CASE/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
[[ -f $LIB1 && -f $LIB2 && -f $FLIB2 && -f $FLIB3 && -f $RLIB ]]
for l in "$LIB1" "$LIB2" "$FLIB2" "$FLIB3" "$RLIB"; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
cp "$CASE/system/fvSolution" fvSolution.orig
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.orig > fvSolution.fixed
echo "fvSolution changes (base -> fixed):"; diff fvSolution.orig fvSolution.fixed || true
[[ $(grep -c 'smoothSolverFixed;' fvSolution.fixed) == 3 && $(grep -c 'PCGFixed;' fvSolution.fixed) == 1 ]]
# grad(p) added to the existing cache block (which already holds grad(U))
awk '{print} /^[[:space:]]*grad\(U\);/ && incache {print "    grad(p);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.fixed > fvSolution.gradp
[[ $(grep -c 'grad(p);' fvSolution.gradp) == 1 ]]
echo "fvSolution change (fixed -> gradp):"; diff fvSolution.fixed fvSolution.gradp || true
awk '{print} /^[[:space:]]*grad\(p\);/ && incache {print "    grad(k);"; print "    grad(omega);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.gradp > fvSolution.kw
[[ $(grep -c 'grad(k);' fvSolution.kw) == 1 && $(grep -c 'grad(omega);' fvSolution.kw) == 1 ]]
echo "fvSolution change (gradp -> kw):"; diff fvSolution.gradp fvSolution.kw || true
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

run() {   # $1 arm (base|cachekw|all|nokomega|nogs)
    local a=$1 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    cd "$CASE"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    local pre=()
    if [[ $a == base || $a == cachekw ]]; then
        if [[ $a == base ]]; then cp "$ROOT/fvSolution.gradp" system/fvSolution; else cp "$ROOT/fvSolution.kw" system/fvSolution; fi
        echo "libs (\"$LIB1\");" >> system/controlDict
        pre=(-genv LD_PRELOAD "$FLIB2" -genv FUSED_VECTOR_SOLVE_REPORT 1)
    else
        cp "$ROOT/fvSolution.kw" system/fvSolution
        echo "libs (\"$LIB2\");" >> system/controlDict
        pre=(-genv LD_PRELOAD "$FLIB3:$RLIB" -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1)
        if [[ $a == nokomega ]]; then pre+=(-genv STEP_REUSE_OFF komega); fi
        if [[ $a == nogs ]]; then pre+=(-genv STEP_REUSE_OFF gs); fi
    fi
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
    echo "ACTIVE arm=$a: $(grep -E 'fusedVectorSolve|stepReuse' "$dir/simpleFoam.log" | sort -u | tr '\n' ' ')"
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
read -r -a ARMS <<< "${ARMS:-base cachekw all}"   # attribution arms: nokomega nogs
[[ ${ARMS[0]} == base ]]
for a in "${ARMS[@]}"; do run "$a"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"   # case restored

for a in "${ARMS[@]:1}"; do
    if diff <(grep -E 'Solving for|continuity errors' run_base/simpleFoam.log) <(grep -E 'Solving for|continuity errors' run_$a/simpleFoam.log) > solver_$a.diff; then
        echo "CHECK_SOLVER_LINES_IDENTICAL $a: yes ($(grep -c 'Solving for' run_base/simpleFoam.log) lines)"
    else
        echo "CHECK_SOLVER_LINES_IDENTICAL $a: NO"; head -6 solver_$a.diff
    fi
    if cmp -s <(grep -v '^#' run_base/postProcessing/forceCoeffsAll/*/coefficient*.dat) <(grep -v '^#' run_$a/postProcessing/forceCoeffsAll/*/coefficient*.dat); then
        echo "CHECK_FORCE_COEFFS_IDENTICAL $a: yes"
    else
        echo "CHECK_FORCE_COEFFS_IDENTICAL $a: NO"
    fi
done
python3 - "$ROOT" "$STEPS" <<'PY'
import re, sys
from pathlib import Path
root, steps = Path(sys.argv[1]), sys.argv[2]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
arms = [a for a in ['base', 'cachekw', 'all', 'nokomega', 'nogs'] if (root / f'run_{a}').is_dir()]
label = {'base': 'base (plugin v1 + fused v2 + grad(p) cache)', 'cachekw': 'base + cache grad(k), grad(omega)',
         'all': 'all (cache + stepReuse + fixedIter v2 + fused v3)',
         'nokomega': 'all without kOmegaSST reuse', 'nogs': 'all without GaussSeidel fold'}
t = {a: float(pat.findall((root / f'run_{a}' / 'simpleFoam.log').read_text())[-1]) for a in arms}
lines = [f'# stepReuse bundle vs current Bonus, 4 nodes x 104, occDrivAerStaticMesh 65M, {steps} steps, same job', '',
         '| Arm | s/step | vs base |', '|---|---:|---:|']
for a in arms:
    lines.append(f"| {label[a]} | {t[a]:.4f} | {'-' if a == 'base' else '%+.2f%%' % (100 * (t[a] / t['base'] - 1))} |")
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n'); print('\n'.join(lines))
PY
echo "STEPREUSE_COMPLETE root=$ROOT"
