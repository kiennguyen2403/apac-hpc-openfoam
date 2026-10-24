#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:25:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N pagg_4n
#PBS -j oe

# GAMG processor agglomeration for p (BONUS F2: allowed fvSolution entry; changes rounding, so Cd check),
# with the current setup (simpleFoamFused + fixedIterSolvers v2 + fusedVectorSolve v3 + production stepReuse
# + cache grad(p), grad(k), grad(omega); measured Bonus 1.2675 s/step, job 182135721).
# APS (job 182137769): per step ~1489 Allreduce and ~1248 halo-exchange rounds, most from the p GAMG coarse
# levels and the 35-iteration coarsest PCG over all 416 ranks.
# Paired test, 4 nodes x 104 ranks, prepared case416_w05, 20 steps per arm, same nodes:
#   base  : no processorAgglomerator
#   mc104 : masterCoarsest, nProcessorsPerMaster 104 (one master per node)
#   mc26  : masterCoarsest, nProcessorsPerMaster 26
#   pf50  : procFaces, nAgglomeratingCells 50
#   pf400 : procFaces, nAgglomeratingCells 400

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
CASE=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
T=/scratch/ya41/kn3896/openfoam/tools
LIB2=$T/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB3=$T/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
RLIB=$T/stepReuse/lib/linux64IcxDPInt32OptSPR/libstepReuse.so
SFF=$T/simpleFoamFused/bin/linux64IcxDPInt32OptSPR/simpleFoamFused

ROOT=/scratch/ya41/kn3896/openfoam/contest-results/pagg/bench_4node_${PBS_JOBID%%.*}
STEPS=20
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "PAGG_FAILED line=$LINENO status=$?"; cp "$ROOT/fvSolution.orig" "$CASE/system/fvSolution" 2>/dev/null || true' ERR
cd "$ROOT"
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/logFiles/$l.log"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
ls "$CASE/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
[[ -f $LIB2 && -f $FLIB3 && -f $RLIB && -x $SFF ]]
for l in "$LIB2" "$FLIB3" "$RLIB" "$SFF"; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
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

# fvSolution per arm: the base one plus the agglomerator entries in solvers.p
for a in mc104 mc26 pf50 pf400; do
    cp "$ROOT/fvSolution.kw" "$ROOT/fvSolution.$a"
    case $a in
        mc104) foamDictionary -entry solvers.p.processorAgglomerator -add masterCoarsest "$ROOT/fvSolution.$a" > /dev/null
               foamDictionary -entry solvers.p.nProcessorsPerMaster -add 104 "$ROOT/fvSolution.$a" > /dev/null ;;
        mc26)  foamDictionary -entry solvers.p.processorAgglomerator -add masterCoarsest "$ROOT/fvSolution.$a" > /dev/null
               foamDictionary -entry solvers.p.nProcessorsPerMaster -add 26 "$ROOT/fvSolution.$a" > /dev/null ;;
        pf50)  foamDictionary -entry solvers.p.processorAgglomerator -add procFaces "$ROOT/fvSolution.$a" > /dev/null
               foamDictionary -entry solvers.p.nAgglomeratingCells -add 50 "$ROOT/fvSolution.$a" > /dev/null ;;
        pf400) foamDictionary -entry solvers.p.processorAgglomerator -add procFaces "$ROOT/fvSolution.$a" > /dev/null
               foamDictionary -entry solvers.p.nAgglomeratingCells -add 400 "$ROOT/fvSolution.$a" > /dev/null ;;
    esac
    echo "fvSolution change (base -> $a):"; diff "$ROOT/fvSolution.kw" "$ROOT/fvSolution.$a" || true
done

run() {   # $1 arm (base|mc104|mc26|pf50|pf400)
    local a=$1 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    cd "$CASE"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    local pre=() app=$SFF lib=$RLIB
    local fvs=$ROOT/fvSolution.kw
    if [[ $a != base ]]; then fvs=$ROOT/fvSolution.$a; fi
    cp "$fvs" system/fvSolution
    echo "libs (\"$LIB2\");" >> system/controlDict
    pre=(-genv LD_PRELOAD "$FLIB3:$lib" -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1 -genv SIMPLEFOAM_FUSED_REPORT 1)
    echo "RUN_START arm=$a utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 "${glibcT[@]}" "${pre[@]}" "$app" -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cp "$ROOT/fvSolution.orig" system/fvSolution
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "ACTIVE arm=$a: $(grep -E 'fusedVectorSolve|stepReuse|simpleFoamFused' "$dir/simpleFoam.log" | sort -u | tr '\n' ' ')"
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
read -r -a ARMS <<< "${ARMS:-base mc104 mc26 pf50 pf400}"
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
arms = [a for a in ['base', 'mc104', 'mc26', 'pf50', 'pf400'] if (root / f'run_{a}').is_dir()]
label = {'base': 'base (no processorAgglomerator)', 'mc104': 'masterCoarsest, 104 per master',
         'mc26': 'masterCoarsest, 26 per master', 'pf50': 'procFaces, 50 cells', 'pf400': 'procFaces, 400 cells'}
t = {a: float(pat.findall((root / f'run_{a}' / 'simpleFoam.log').read_text())[-1]) for a in arms}
def cd_of(a):
    f = sorted((root / f'run_{a}' / 'postProcessing' / 'forceCoeffsAll').glob('*/coefficient*.dat'))[0]
    last = [l for l in f.read_text().splitlines() if l.strip() and not l.startswith('#')][-1]
    return float(last.split()[1])
cd = {a: cd_of(a) for a in arms}
lines = [f'# GAMG processor agglomeration for p, current Bonus setup, 4 nodes x 104, occDrivAerStaticMesh 65M, {steps} steps, same job', '',
         '| Arm | s/step | vs base | Cd step %s | Cd vs base |' % steps, '|---|---:|---:|---:|---:|']
for a in arms:
    lines.append(f"| {label[a]} | {t[a]:.4f} | {'-' if a == 'base' else '%+.2f%%' % (100 * (t[a] / t['base'] - 1))} | {cd[a]:.8f} | {'-' if a == 'base' else '%+.4f%%' % (100 * (cd[a] / cd['base'] - 1))} |")
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n'); print('\n'.join(lines))
PY
echo "PAGG_COMPLETE root=$ROOT"
