#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=208,mem=400gb
#PBS -l walltime=00:20:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N poll_2n
#PBS -j oe

# nPollProcInterfaces=-1 vs default on 2 nodes x 104 ranks (inter-node neighbours over InfiniBand),
# occDrivAerStaticMesh 65M, 40 steps. Switch set with -opt-switch (no case file edited for it).
# Reuses the prepared 208-rank case of the IPM pilot (plain scotch, uncollated) -- same case for
# both runs, so only the switch differs. ICX OptSPR + Intel MPI, team scatter pinning, 2-node
# tuned Allreduce file. One run each: base, then poll.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
NP=208
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
TEAM=/scratch/ya41/apac-2026/openfoam
P=/scratch/ya41/kn3896/openfoam/contest-results/communication/pilot_ipm_contest_2n_15iter_178472991.gadi-pbs
C=$P/case
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/optswitch/poll_2node_${PBS_JOBID%%.*}
STEPS=40

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "POLL2N_FAILED line=$LINENO status=$?"; sed -i "s/^endTime .*/endTime         15;/" "$C/system/controlDict"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 2 ]]
cp "$TEAM/run/2node/intel-scatter-cpus.awk" "$TEAM/run/2node/intelmpi-allreduce-2n-104ppn-scatter.dat" .
for l in decomposePar renumberMesh potentialFoam applyBoundaryLayer; do grep -q '^End' "$P/preparation/reused-source/$l.log"; done
[[ $(find "$C" -maxdepth 1 -type d -name 'processor[0-9]*' | wc -l) == $NP ]]
cmp "$C/system/fvSolution" "$C/system/fvSolution.fixedIter"
ls "$C/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
if [[ -d $C/postProcessing ]]; then mv "$C/postProcessing" "$ROOT/postProcessing.before"; fi
cp "$C/system/controlDict" "$ROOT/controlDict.original"
sed -i "s/^endTime .*/endTime         $STEPS;/" "$C/system/controlDict"
[[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' "$C/system/controlDict") == $STEPS ]]

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
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-2n-104ppn-scatter.dat
mpirun -np $NP -ppn 104 hostname | sort | uniq -c > placement.txt
if [[ $(wc -l < placement.txt) != 2 ]] || awk '$1 != 104 {bad=1} END {exit !bad}' placement.txt; then
    echo "PLACEMENT_WRONG:"; cat placement.txt; exit 1
fi

run() {   # $1 variant (base|poll), $2 repeat
    local v=$1 r=$2 dir=$ROOT/run_${1}_$2 status=0 sw=()
    [[ $v == poll ]] && sw=(-opt-switch nPollProcInterfaces=-1)
    mkdir "$dir"
    cd "$C"
    echo "RUN_START variant=$v repeat=$r utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np $NP -ppn 104 simpleFoam -parallel "${sw[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS variant=$v repeat=$r Cd_step$STEPS=$cd polling=$(grep -m1 -oE 'polling iterations : -?[0-9]+' "$dir/simpleFoam.log" | awk '{print $NF}') $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

run base 1
run poll 1
cp "$ROOT/controlDict.original" "$C/system/controlDict"   # leave the pilot case as found

python3 - "$ROOT" <<'EOF'
import re, sys
from pathlib import Path
root = Path(sys.argv[1])
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
t = lambda v, r: float(pat.findall((root / f'run_{v}_{r}' / 'simpleFoam.log').read_text())[-1])
m = {v: t(v, 1) for v in ('base', 'poll')}
lines = ['# nPollProcInterfaces=-1 vs default, 2 nodes x 104, occDrivAerStaticMesh 65M, 40 steps', '',
         'ICX OptSPR + Intel MPI, scatter pinning, 2-node tuned Allreduce; IPM-pilot 208-rank scotch case. One run each, same job.', '',
         '| Variant | s/step | vs base |', '|---|---:|---:|']
for v in ('base', 'poll'):
    rel = f'{100*(m[v]/m["base"]-1):+.2f}%' if v != 'base' else '-'
    lines.append(f'| {v} | {m[v]:.4f} | {rel} |')
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "POLL2N_COMPLETE root=$ROOT"
