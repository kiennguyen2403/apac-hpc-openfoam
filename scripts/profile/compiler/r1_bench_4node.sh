#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N r1_bench_4n
#PBS -j oe

# Decomposition search round 1, step 3: paired 4-node test of the top-3 scored decompositions
# against the current best, 10 steps each, same job / nodes, one run each (no mirror).
#   cur     : case416_w05, scotch (4 104)/(1 0.5)  (prepared, reused)
#   sock211 : (4 2 52)/(2 1 1)       predicted -2.57%
#   w5      : (4 104)/(5 1)          predicted -2.29%
#   (sbal01 decomposed but dropped from the run to save SU)
# Parallel prep of the three new cases with the team chain (site OpenFOAM/2512): renumberMesh,
# potentialFoam, applyBoundaryLayer; reused if complete. Solver: current best setup = ICX OptSPR +
# Intel MPI, scatter pinning, 4-node tuned Allreduce, glibcT. Accept rule: >= 0.7% faster than cur.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
TIMER_SRC=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case/dynamicCode
SEARCH=/scratch/ya41/kn3896/openfoam/contest-results/decomp-search
ROOT=$SEARCH/r1_bench_${PBS_JOBID%%.*}
STEPS=10
ARMS=(cur sock211 w5)
declare -A CASE=([cur]=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
                 [sock211]=$SEARCH/case416_sock211 [w5]=$SEARCH/case416_w5 [sbal01]=$SEARCH/case416_sbal01)
declare -A PRED=([cur]=0 [sock211]=-2.57 [w5]=-2.29 [sbal01]=-1.67)
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "R1_BENCH_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for a in "${ARMS[@]}"; do grep -q '^End' "${CASE[$a]}/logFiles/10_decomposePar.log"; done

# --- parallel prep of new cases (team chain, site module) ------------------------------
for a in sock211 w5; do
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
done

# --- solver environment (current best setup) ---------------------------------------
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
[[ $(tr ',' '\n' <<< "$I_MPI_PIN_PROCESSOR_LIST" | wc -l) == 104 ]]
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
    [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == $STEPS ]]
    echo "RUN_START arm=$a utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 600s mpirun -np 416 -ppn 104 "${glibcT[@]}" simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    cd "$ROOT"
    if [[ $status != 0 ]] || ! grep -qx End "$dir/simpleFoam.log" || ! grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"; then
        echo "RUN_FAILED arm=$a status=$status: $(grep -m1 -iE 'FATAL|error|signal' "$dir/simpleFoam.log" | cut -c1-160)"
        return 0
    fi
    echo "RUN_SUCCESS arm=$a $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
for a in "${ARMS[@]}"; do run "$a"; done

module load python3/3.12.1 > /dev/null 2>&1
python3 -I - "$ROOT" "$STEPS" "${ARMS[*]}" "$(for a in "${ARMS[@]}"; do printf '%s ' "${PRED[$a]}"; done)" <<'PY'
import re, sys, statistics as st
from pathlib import Path
root, steps, arms, preds = Path(sys.argv[1]), int(sys.argv[2]), sys.argv[3].split(), [float(x) for x in sys.argv[4].split()]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
res = {}
for a in arms:
    log = root / f'run_{a}' / 'simpleFoam.log'
    txt = log.read_text() if log.exists() else ''
    m = pat.findall(txt)
    ex = [float(x) for x in re.findall(r'^ExecutionTime = ([0-9.]+) s', txt, re.M)]
    d = [b - a_ for a_, b in zip(ex, ex[1:])]
    sd = 100 * st.stdev(d) / st.mean(d) if len(d) > 2 else float('nan')
    rows = [l.split() for p in (root / f'run_{a}' / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    res[a] = (float(m[-1]) if m else None, float(rows[-1][1]) if rows else None, sd)
base = res['cur'][0]
out = [f'# Decomposition search round 1: paired 4-node test, {steps} steps each, same job', '',
       'ICX OptSPR + Intel MPI, scatter pinning, 4-node tuned Allreduce, glibcT. One run each.',
       'Accept: >= 0.7% faster than cur (paired noise at 10-15 steps ~0.35%).', '',
       f'| Arm | predicted | s/step | measured vs cur | per-step sd | Cd step {steps} |', '|---|---:|---:|---:|---:|---:|']
for a, p in zip(arms, preds):
    t, cd, sd = res[a]
    rel = f'{100 * (t / base - 1):+.2f}%' if (t and base) else 'n/a'
    out.append(f"| {a} | {p:+.2f}% | {t:.4f} | {rel} | {sd:.2f}% | {cd:.5f} |" if t else f'| {a} | {p:+.2f}% | failed | | | |')
(root / 'REPORT.md').write_text('\n'.join(out) + '\n')
print('\n'.join(out))
PY
echo "R1_BENCH_COMPLETE root=$ROOT"
