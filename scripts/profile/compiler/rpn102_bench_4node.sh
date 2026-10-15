#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N rpn102_4n
#PBS -j oe

# 102 vs 104 ranks/node on 4 nodes, occDrivAerStaticMesh 65M, 40-step screens.
#   arm 416 (104/node): an0642's prepared scotch domains (4 104) case, used read-only via links
#   arm 408 (102/node): case408 from rpn102_decompose.sh; parallel prep done here (site module,
#                       as the team chain) -- renumberMesh, potentialFoam, applyBoundaryLayer
# Solver: custom ICX OptSPR + Intel MPI 2021.17.2, team scatter pinning (first P entries of the
# scatter list) and the 4-node tuned Allreduce file. Mirrored order 416, 408, 408, 416.
# The compiled timer is copied in, so no warmup runs are needed.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
A416=/scratch/ya41/an0642/case/scotch_4_node_multi
TIMER_SRC=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case/dynamicCode
BASE=/scratch/ya41/kn3896/openfoam/contest-results/ranks-per-node
ROOT=$BASE/bench_${PBS_JOBID%%.*}
C408=$BASE/case408
C416=$ROOT/case416
fh416=(-fileHandler hostCollated)
fh408=(-fileHandler collated -ioRanks '(0 102 204 306)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "RPN_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .

# --- arm 416: read-only links to an0642's prepared case --------------------------
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$A416/logFiles/$l.log"; done
mkdir -p "$C416/constant"
cp -a "$A416/system" "$C416/"
cp "$TEAM/dicts/4node/controlDict" "$C416/system/controlDict"     # team file, endTime 40
cp -a "$A416"/constant/{transportProperties,turbulenceProperties} "$C416/constant/"
ln -s "$MESH" "$C416/constant/polyMesh"
for p in "$A416"/processors416*; do ln -s "$p" "$C416/${p##*/}"; done
cp -a "$TIMER_SRC" "$C416/dynamicCode"

# --- arm 408: parallel prep (site module, like the team chain) ---------------------
grep -q '^End' "$C408/logFiles/10_decomposePar.log"
prep_done=1
for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$C408/logFiles/$l.log" 2>/dev/null || prep_done=0; done
if [[ $prep_done == 1 ]]; then echo "REUSING prepared 408 case $C408"; else
(
    set +eu
    module purge
    module load OpenFOAM/2512
    set -eu
    cd "$C408"
    par=(mpirun -np 408 --map-by ppr:102:node --bind-to core)
    "${par[@]}" renumberMesh -constant -overwrite -parallel "${fh408[@]}" > logFiles/20_renumberMesh.log 2>&1
    "${par[@]}" potentialFoam -initialiseUBCs -parallel "${fh408[@]}" > logFiles/30_potentialFoam.log 2>&1
    "${par[@]}" applyBoundaryLayer -ybl "0.0450244" -parallel "${fh408[@]}" > logFiles/40_applyBoundaryLayer.log 2>&1
    for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' logFiles/$l.log; done
)
fi
[[ -d $C408/dynamicCode ]] || cp -a "$TIMER_SRC" "$C408/dynamicCode"
echo "PREPARED_408 $(grep -E 'Max number of processor patches' "$C408/logFiles/10_decomposePar.log")"
echo "REFERENCE_416 $(grep -E 'Max number of processor patches' "$A416/logFiles/10_decomposePar.log")"

# --- solver environment ----------------------------------------------------------
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
# Under PBS, Intel MPI follows the scheduler's placement and ignores -ppn by default: job
# 180765170 got 104/104/104/96 ranks for -ppn 102 (two ranks per core, ~50x slower).
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-4n-104ppn-scatter.dat
SCATTER=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
[[ $(tr ',' '\n' <<< "$SCATTER" | wc -l) == 104 ]]
echo "simpleFoam=$(command -v simpleFoam)" > build-metadata.txt

run() {   # $1 ranks (416|408), $2 label
    local np=$1 label=$2 ppn cdir dir
    if [[ $np == 416 ]]; then ppn=104 cdir=$C416; else ppn=102 cdir=$C408; fi
    dir=$ROOT/run_${np}_$label
    mkdir "$dir"
    local fh=("${fh416[@]}"); [[ $np == 408 ]] && fh=("${fh408[@]}")
    export I_MPI_PIN_PROCESSOR_LIST=$(tr ',' '\n' <<< "$SCATTER" | head -n $ppn | paste -sd,)
    echo "$I_MPI_PIN_PROCESSOR_LIST" > "$dir/pin-list.txt"
    cd "$cdir"
    [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == 40 ]]
    mpirun -np $np -ppn $ppn hostname | sort | uniq -c > "$dir/placement.txt"
    if [[ $(wc -l < "$dir/placement.txt") != 4 ]] || awk -v p=$ppn '$1 != p {bad=1} END {exit !bad}' "$dir/placement.txt"; then
        echo "PLACEMENT_WRONG ranks=$np ppn=$ppn:"; cat "$dir/placement.txt"; return 1
    fi
    echo "RUN_START ranks=$np ppn=$ppn run=$label utc=$(date -u +%FT%TZ)"
    local status=0
    timeout --kill-after=30s 900s mpirun -np $np -ppn $ppn simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == 40 ]]
    grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS ranks=$np run=$label Cd_step40=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

# 408 arm only; the 416 reference is run_416_t1 of job 180765170 (same nodes type, same build/settings).
run 408 t1
run 408 t2
REF416=$BASE/bench_180765170/run_416_t1
mkdir -p "$ROOT/run_416_t1" && cp "$REF416/simpleFoam.log" "$ROOT/run_416_t1/"

python3 - "$ROOT" <<'EOF'
import re, sys
from pathlib import Path
root = Path(sys.argv[1])
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
t = lambda n, r: float(pat.findall((root / f'run_{n}_{r}' / 'simpleFoam.log').read_text())[-1])
m = {416: (t(416, 't1'), t(416, 't1'), float('nan')),
     408: ((t(408, 't1') + t(408, 't2')) / 2, t(408, 't1'), t(408, 't2'))}
lines = ['# 102 vs 104 ranks/node, 4 nodes, occDrivAerStaticMesh 65M, 40 steps', '',
         'ICX OptSPR + Intel MPI 2021.17.2, scatter pinning, 4-node tuned Allreduce file.',
         '416 = single run from job 180765170 (reference); 408 = two runs in this job.', '',
         '| Ranks | Per node | Repeat 1 (s) | Repeat 2 (s) | Mean (s/step) | vs 416 |', '|---:|---:|---:|---:|---:|---:|']
for n, ppn in ((416, 104), (408, 102)):
    mean, a, b = m[n]
    lines.append(f'| {n} | {ppn} | {a:.4f} | {b:.4f} | {mean:.4f} | {100 * (mean / m[416][0] - 1):+.2f}% |')
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "RPN_COMPLETE root=$ROOT"
