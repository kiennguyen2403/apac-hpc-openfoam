#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N scalecorr150
#PBS -j oe

# BONUS (fvSolution numerics change, not bit-identical): GAMG scaleCorrection off for p. With
# scaling on (default for symmetric matrices) every V-cycle level does an extra A*x pass, a global
# reduction and a Jacobi-like update of the correction. Full-length check: 4 nodes x 104, 150 steps,
# noscale arm only, Cd compared with the Basic/Bonus Cd history of job 181332282 (run_basic_1). Arm = current Bonus setup (fixedIterSolvers plugin + fused vector solve via
# LD_PRELOAD, glibcT, ICX OptSPR + Intel MPI, scatter pinning, tuned Allreduce). Runs in a linked
# copy of case416_w05 (read-only processor directories) so the queued final-runs job is unaffected.
#   bonus    : p GAMG as in the team fvSolution
#   noscale  : same + scaleCorrection false;

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
CUR=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/scalecorr/full150_4node_${PBS_JOBID%%.*}
REFRUN=/scratch/ya41/kn3896/openfoam/contest-results/final/runs_4n_181332282/run_basic_1
CASEDIR=$ROOT/case
STEPS=150
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "SCALECORR_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CUR/logFiles/$l.log"; done
[[ -f $LIB && -f $FLIB ]]
cmp "$CUR/system/fvSolution.fixedIter" /scratch/ya41/apac-2026/openfoam/run/case/occDrivAerStaticMesh/system/fvSolution.fixedIter
cp "$CUR/system/fvSolution.fixedIter" fvSolution.basic
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.basic > fvSolution.bonus
[[ $(grep -c 'smoothSolverFixed;' fvSolution.bonus) == 3 && $(grep -c 'PCGFixed;' fvSolution.bonus) == 1 ]]
# scaleCorrection false in the p solver only (the line after "solver GAMG;")
awk '{print} /^[[:space:]]*solver[[:space:]]+GAMG;/ {print "        scaleCorrection false;"}' fvSolution.bonus > fvSolution.noscale
[[ $(grep -c 'scaleCorrection false;' fvSolution.noscale) == 1 ]]
echo "fvSolution change (bonus -> noscale):"; diff fvSolution.bonus fvSolution.noscale || true

# Both arms get their own case directory (case416_w05 itself is in use by the final-runs job):
# our case setup + links to prepared, read-only processor directories
mklinked() {   # $1 new case dir, $2 source of processor directories
    mkdir -p "$1/constant"
    cp -a "$CUR/system" "$1/"
    cp "$ROOT/fvSolution.basic" "$1/system/fvSolution"
    cp -a "$CUR"/constant/{transportProperties,turbulenceProperties} "$1/constant/"
    ln -s "$MESH" "$1/constant/polyMesh"
    for p in "$2"/processors416_*; do ln -s "$p" "$1/${p##*/}"; done
    cp -a "$CUR/dynamicCode" "$1/dynamicCode"
    echo "$1: $(ls -d "$1"/processors416_* | wc -l) processor dirs linked"
}
mklinked "$CASEDIR" "$CUR"

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
mpirun -np 416 -ppn 104 hostname | sort | uniq -c > placement.txt
if [[ $(wc -l < placement.txt) != 4 ]] || awk '$1 != 104 {bad=1} END {exit !bad}' placement.txt; then
    echo "PLACEMENT_WRONG:"; cat placement.txt; exit 1
fi

run() {   # $1 arm (bonus|noscale)
    local a=$1 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    cd "$CASEDIR"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    echo "libs (\"$LIB\");" >> system/controlDict
    cp "$ROOT/fvSolution.$a" system/fvSolution
    echo "RUN_START arm=$a utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 "${glibcT[@]}" -genv LD_PRELOAD "$FLIB" -genv FUSED_VECTOR_SOLVE_REPORT 1 simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cp "$ROOT/fvSolution.basic" system/fvSolution
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    grep -q "fusedVectorSolve: U fused path" "$dir/simpleFoam.log"
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
run noscale

python3 - "$ROOT" "$REFRUN" <<'PY'
import re, sys
from pathlib import Path
root, ref = Path(sys.argv[1]), Path(sys.argv[2])
def series(d):
    rows = [l.split() for p in (d / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return {int(float(r[0])): float(r[1]) for r in rows}
n, b = series(root / 'run_noscale'), series(ref)
t = float(re.findall(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)', (root / 'run_noscale' / 'simpleFoam.log').read_text())[-1])
out = ['# GAMG scaleCorrection off for p: 150-step Cd check, 4 nodes x 104, Bonus setup', '',
       f'noscale s/step = {t:.4f} (Bonus with fused solve, job 181334343: 1.6066; Basic, job 181332282: 1.9765)', '',
       '| Step | Cd Basic/Bonus (job 181332282) | Cd noscale | difference |', '|---:|---:|---:|---:|']
for st in (10, 20, 40, 80, 120, 150):
    if st in n and st in b:
        out.append(f'| {st} | {b[st]:.8f} | {n[st]:.8f} | {100 * (n[st] / b[st] - 1):+.3f}% |')
final = 100 * (n[150] / b[150] - 1)
out += ['', f'Bonus Cd rule (within +-0.5% of baseline at step 150): **{"PASS" if abs(final) <= 0.5 else "FAIL"}** ({final:+.3f}%)']
(root / 'REPORT.md').write_text('\n'.join(out) + '\n'); print('\n'.join(out))
PY
echo "SCALECORR_COMPLETE root=$ROOT"
