#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:50:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N weighted_4n
#PBS -j oe

# BONUS: weighted scotch decompositions (weighted_decompose_416.sh) vs the current scotch
# (4 104)/(1 0.5), paired, 4 nodes x 104, 20 steps each, same job, Bonus setup (fixedIterSolvers
# plugin + fused vector solve v2, glibcT, ICX OptSPR + Intel MPI, scatter pinning, tuned
# Allreduce). Parallel prep of the weighted cases with the team chain (site OpenFOAM/2512).
#   cur   : case416_w05
#   wface : weights from cell faces (4*internal + 7*wall)
#   wmeas : faces weights scaled per rank by measured cost (diagnosis job 181350543)

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
TIMER_SRC=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case/dynamicCode
W=/scratch/ya41/kn3896/openfoam/contest-results/weighted
ROOT=$W/bench_4node_${PBS_JOBID%%.*}
STEPS=20
ARMS=(cur wface wmeas)
declare -A CASE=([cur]=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
                 [wface]=$W/case416_wface [wmeas]=$W/case416_wmeas)
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "WEIGHTED_BENCH_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for a in "${ARMS[@]}"; do grep -q '^End' "${CASE[$a]}/logFiles/10_decomposePar.log"; done
[[ -f $LIB && -f $FLIB ]]
cp "${CASE[cur]}/system/fvSolution.fixedIter" fvSolution.basic
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.basic > fvSolution.bonus
[[ $(grep -c 'smoothSolverFixed;' fvSolution.bonus) == 3 && $(grep -c 'PCGFixed;' fvSolution.bonus) == 1 ]]

# --- parallel prep of the weighted cases (team chain, site module) ---------------------
for a in wface wmeas; do
    C=${CASE[$a]}
    done_=1
    for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$C/logFiles/$l.log" 2>/dev/null || done_=0; done
    if [[ $done_ == 1 ]]; then echo "REUSING prepared $a"; continue; fi
    (
        set +eu
        module purge
        module load OpenFOAM/2512
        set -eu
        cd "$C"
        par=(mpirun -np 416 --map-by ppr:104:node --bind-to core)
        echo "PREP_START $a utc=$(date -u +%FT%TZ)"
        "${par[@]}" renumberMesh -constant -overwrite -parallel "${fh[@]}" > logFiles/20_renumberMesh.log 2>&1
        "${par[@]}" potentialFoam -initialiseUBCs -parallel "${fh[@]}" > logFiles/30_potentialFoam.log 2>&1
        "${par[@]}" applyBoundaryLayer -ybl "0.0450244" -parallel "${fh[@]}" > logFiles/40_applyBoundaryLayer.log 2>&1
        for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' logFiles/$l.log; done
    )
done
for a in "${ARMS[@]}"; do
    [[ -d ${CASE[$a]}/dynamicCode ]] || cp -a "$TIMER_SRC" "${CASE[$a]}/dynamicCode"
    ls "${CASE[$a]}/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
    cmp "${CASE[$a]}/system/fvSolution" "${CASE[$a]}/system/fvSolution.fixedIter"
    echo "$a: $(grep -E 'Max number of (cells|processor patches|faces between)' "${CASE[$a]}/logFiles/10_decomposePar.log" | tr -s ' ' | tr '\n' ' ')"
done

# --- solver environment (Bonus setup) --------------------------------------------------
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

run() {   # $1 arm
    local a=$1 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    cd "${CASE[$a]}"
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/postProcessing.before"; fi
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    echo "libs (\"$LIB\");" >> system/controlDict
    cp "$ROOT/fvSolution.bonus" system/fvSolution
    echo "RUN_START arm=$a utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 "${glibcT[@]}" -genv LD_PRELOAD "$FLIB" -genv FUSED_VECTOR_SOLVE_REPORT 1 \
        simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cp "$ROOT/fvSolution.basic" system/fvSolution
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    grep -q 'fusedVectorSolve: U fused path' "$dir/simpleFoam.log"
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
for a in "${ARMS[@]}"; do run "$a"; done

python3 - "$ROOT" "$STEPS" <<'PY'
import re, sys
from pathlib import Path
root, steps = Path(sys.argv[1]), sys.argv[2]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
def t(a): return float(pat.findall((root / f'run_{a}' / 'simpleFoam.log').read_text())[-1])
def cd(a):
    rows = [l.split() for p in (root / f'run_{a}' / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return float(rows[-1][1])
c = t('cur')
out = [f'# Weighted scotch decompositions vs current, 4 nodes x 104, {steps} steps, Bonus setup, same job', '',
       '| Arm | s/step | vs cur | Cd step ' + steps + ' | Cd vs cur |', '|---|---:|---:|---:|---:|']
cd0 = cd('cur')
for a, label in (('cur', 'scotch (4 104)/(1 0.5), unweighted'), ('wface', 'weighted: faces'), ('wmeas', 'weighted: measured')):
    x, y = t(a), cd(a)
    dt = '-' if a == 'cur' else '%+.2f%%' % (100 * (x / c - 1))
    dc = '-' if a == 'cur' else '%+.4f%%' % (100 * (y / cd0 - 1))
    out.append('| %s | %.4f | %s | %.8f | %s |' % (label, x, dt, y, dc))
(root / 'REPORT.md').write_text('\n'.join(out) + '\n'); print('\n'.join(out))
PY
echo "WEIGHTED_BENCH_COMPLETE root=$ROOT"
