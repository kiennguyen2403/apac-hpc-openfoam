#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N numaflat_4n
#PBS -j oe

# Topology-aware rank mapping vs current, 4 nodes x 104 ranks, occDrivAerStaticMesh 65M, 40 steps.
#   cur      : scotch domains (4 104), weights (1 0.5) + SCATTER pinning (team setting)
#              -- prepared case416_w05, reused
#   numaflat : scotch domains (4 2 4 13), weights (1 0.5 0.5 0.5) + COMPACT pinning
#              (rank i on core i: NUMA k = cores 13k..13k+12, sockets 0-51 / 52-103)
#              -- case416_numaflat from numaflat_decompose.sh; parallel prep here (site module)
# Same job/nodes, ICX OptSPR + Intel MPI, 4-node tuned Allreduce file, glibc. Basic-legal:
# scotch coefficients + runtime pinning only. Cd must stay within 0.43-0.44 at step 40.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
BASE=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen
CCUR=$BASE/case416_w05
CNF=$BASE/case416_numaflat
TIMER_SRC=$CCUR/dynamicCode
ROOT=$BASE/bench_numaflat_${PBS_JOBID%%.*}
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')
STEPS=40

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "NUMAFLAT_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CCUR/logFiles/$l.log"; done
grep -q '^nDomainWeights (1 0.5);' "$CCUR/system/include/caseDefinition"
grep -q '^End' "$CNF/logFiles/10_decomposePar.log"
grep -q '^nDomains (4 2 4 13);' "$CNF/system/include/caseDefinition"
grep -q '^nDomainWeights (1 0.5 0.5 0.5);' "$CNF/system/include/caseDefinition"
for c in "$CCUR" "$CNF"; do cmp "$c/system/fvSolution" "$c/system/fvSolution.fixedIter"; done

# --- numaflat: parallel prep (site module, like the team chain), reused if complete -------
prep_done=1
for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CNF/logFiles/$l.log" 2>/dev/null || prep_done=0; done
if [[ $prep_done == 1 ]]; then echo "REUSING prepared numaflat case $CNF"; else
(
    set +eu
    module purge
    module load OpenFOAM/2512
    set -eu
    cd "$CNF"
    par=(mpirun -np 416 --map-by ppr:104:node --bind-to core)
    "${par[@]}" renumberMesh -constant -overwrite -parallel "${fh[@]}" > logFiles/20_renumberMesh.log 2>&1
    "${par[@]}" potentialFoam -initialiseUBCs -parallel "${fh[@]}" > logFiles/30_potentialFoam.log 2>&1
    "${par[@]}" applyBoundaryLayer -ybl "0.0450244" -parallel "${fh[@]}" > logFiles/40_applyBoundaryLayer.log 2>&1
    for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' logFiles/$l.log; done
)
fi
[[ -d $CNF/dynamicCode ]] || cp -a "$TIMER_SRC" "$CNF/dynamicCode"
ls "$CNF/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
echo "PREPARED_NUMAFLAT $(grep -E 'Max number of (cells|processor patches|faces between)' "$CNF/logFiles/10_decomposePar.log" | tr '\n' ' ')"

# --- solver environment (team runtime settings) -------------------------------------------
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
export I_MPI_PIN=1 I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
SCATTER=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
COMPACT=$(seq -s, 0 103)
[[ $(tr ',' '\n' <<< "$SCATTER" | wc -l) == 104 && $(tr ',' '\n' <<< "$COMPACT" | wc -l) == 104 ]]
lscpu > lscpu.txt
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-4n-104ppn-scatter.dat
mpirun -np 416 -ppn 104 hostname | sort | uniq -c > placement.txt
if [[ $(wc -l < placement.txt) != 4 ]] || awk '$1 != 104 {bad=1} END {exit !bad}' placement.txt; then
    echo "PLACEMENT_WRONG:"; cat placement.txt; exit 1
fi

run() {   # $1 arm (cur|numaflat)
    local arm=$1 cdir dir
    dir=$ROOT/run_$arm
    mkdir "$dir"
    if [[ $arm == cur ]]; then
        cdir=$CCUR; export I_MPI_PIN_ORDER=scatter I_MPI_PIN_PROCESSOR_LIST=$SCATTER
    else
        cdir=$CNF; export I_MPI_PIN_ORDER=compact I_MPI_PIN_PROCESSOR_LIST=$COMPACT
    fi
    echo "$I_MPI_PIN_PROCESSOR_LIST" > "$dir/pin-list.txt"
    cd "$cdir"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == $STEPS ]]
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/postProcessing.before"; fi
    echo "RUN_START arm=$arm pin=$(cut -c1-40 "$dir/pin-list.txt")... utc=$(date -u +%FT%TZ)"
    local status=0
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS arm=$arm Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

run cur
run numaflat

python3 - "$ROOT" <<'EOF'
import re, sys
from pathlib import Path
root = Path(sys.argv[1])
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
def t(a): return float(pat.findall((root / f'run_{a}' / 'simpleFoam.log').read_text())[-1])
def cd(a):
    rows = [l.split() for p in (root / f'run_{a}' / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return float(rows[-1][1])
cur, nf = t('cur'), t('numaflat')
lines = ['# Topology-aware rank mapping, 4 nodes x 104, occDrivAerStaticMesh 65M, 40 steps', '',
         'ICX OptSPR + Intel MPI, 4-node tuned Allreduce, glibc. One run each, same job/nodes.', '',
         '| Arm | scotch domains / weights | Pinning | s/step | vs cur | Cd step 40 (0.43-0.44) |',
         '|---|---|---|---:|---:|---:|',
         f'| cur | (4 104) / (1 0.5) | scatter | {cur:.4f} | - | {cd("cur"):.5f} |',
         f'| numaflat | (4 2 4 13) / (1 0.5 0.5 0.5) | compact | {nf:.4f} | {100*(nf/cur-1):+.2f}% | {cd("numaflat"):.5f} |']
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "NUMAFLAT_COMPLETE root=$ROOT"
