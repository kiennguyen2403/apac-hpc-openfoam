#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N optswitch_1n
#PBS -j oe

# OpenFOAM OptimisationSwitches (Basic-allowed runtime switches, rules F1), set on the command line
# with -opt-switch (no case file edited). 1 node, 104 ranks, occDrivAerStaticMesh 65M, 40 steps.
#   base    : defaults
#   poll    : nPollProcInterfaces=-1   (dispatch each processor interface when its data arrives)
#   unified : unifiedGeometricField=1  (fields with flattened boundary storage)
#   both    : both switches
# Reuses the prepared 1-node scotch case (no decomposition). ICX OptSPR + Intel MPI, team scatter
# pinning + 1-node tuned Allreduce. A failed variant is recorded and the job continues.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
NP=104
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
TEAM=/scratch/ya41/apac-2026/openfoam
CASE=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/optswitch/drivaer_1node_${PBS_JOBID%%.*}
fh="-fileHandler hostCollated"
STEPS=40

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "OPTSW_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
cp "$TEAM/run/1node/intel-scatter-cpus.awk" "$TEAM/run/1node/intelmpi-allreduce-1n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/../logFiles/$l.log"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
ls "$CASE/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
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
export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
export I_MPI_PIN_PROCESSOR_LIST
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-1n-104ppn-scatter.dat

switches() {
    case $1 in
        base)    ;;
        poll)    echo "-opt-switch nPollProcInterfaces=-1" ;;
        unified) echo "-opt-switch unifiedGeometricField=1" ;;
        both)    echo "-opt-switch nPollProcInterfaces=-1 -opt-switch unifiedGeometricField=1" ;;
    esac
}

run() {   # $1 variant
    local v=$1 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    cd "$CASE"
    cp "$TEAM/dicts/1node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == $STEPS ]]
    echo "RUN_START variant=$v args='$(switches "$v")' utc=$(date -u +%FT%TZ)"
    # shellcheck disable=SC2046
    timeout --kill-after=30s 1200s mpirun -np $NP simpleFoam -parallel $fh $(switches "$v") > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    if [[ $status != 0 ]] || ! grep -qx End "$dir/simpleFoam.log" || ! grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"; then
        echo "RUN_FAILED variant=$v status=$status: $(grep -m1 -iE 'FATAL|error|signal' "$dir/simpleFoam.log" | cut -c1-160)"
        return 0
    fi
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS variant=$v Cd_step$STEPS=$cd polling=$(grep -m1 -oE 'polling iterations : -?[0-9]+' "$dir/simpleFoam.log" | awk '{print $NF}') $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

for v in base poll unified both; do run "$v"; done
cp "$TEAM/dicts/1node/controlDict" "$CASE/system/controlDict"

python3 - "$ROOT" <<'EOF'
import re, sys
from pathlib import Path
root = Path(sys.argv[1])
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
def t(v):
    f = root / f'run_{v}' / 'simpleFoam.log'
    m = pat.findall(f.read_text()) if f.exists() else []
    return float(m[-1]) if m else None
def cd(v):
    rows = [l.split() for p in (root / f'run_{v}' / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return float(rows[-1][1]) if rows else None
b = t('base')
lines = ['# OptimisationSwitches (-opt-switch), 1 node x 104, occDrivAerStaticMesh 65M, 40 steps', '',
         '| Variant | Switches | s/step | vs base | Cd step 40 (0.43-0.44) |', '|---|---|---:|---:|---:|']
desc = {'base': 'defaults', 'poll': 'nPollProcInterfaces=-1', 'unified': 'unifiedGeometricField=1', 'both': 'both'}
for v in ('base', 'poll', 'unified', 'both'):
    tv, cv = t(v), cd(v)
    if tv is None: lines.append(f'| {v} | {desc[v]} | failed | | |'); continue
    rel = f'{100*(tv/b-1):+.2f}%' if b and v != 'base' else '-'
    lines.append(f'| {v} | {desc[v]} | {tv:.4f} | {rel} | {cv:.5f} |')
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "OPTSW_COMPLETE root=$ROOT"
