#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N dw05_4n
#PBS -j oe

# scotch domainWeights (1 0.5) vs the current (1 0.01), 4 nodes x 104 ranks, occDrivAerStaticMesh
# 65M, 40 steps, one run each in the same job (same nodes). Standard glibc allocator.
#   cur : an0642's prepared scotch domains (4 104) weights (1 0.01) case, read-only via links
#   w05 : case416_w05 from dw05_decompose.sh; parallel prep done here (site module, as the
#         team chain) -- renumberMesh, potentialFoam, applyBoundaryLayer; reused if complete
# Solver: ICX OptSPR + Intel MPI 2021.17.2, team scatter pinning, 4-node tuned Allreduce file.
# The compiled timer is copied in, so no warmup runs are needed.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
A416=/scratch/ya41/an0642/case/scotch_4_node_multi
TIMER_SRC=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case/dynamicCode
BASE=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen
ROOT=$BASE/bench_dw05_${PBS_JOBID%%.*}
CW05=$BASE/case416_w05
CCUR=$ROOT/case416_cur
fhcur=(-fileHandler hostCollated)
fhw05=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "DW05_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .

# --- cur: read-only links to an0642's prepared case --------------------------------
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$A416/logFiles/$l.log"; done
mkdir -p "$CCUR/constant"
cp -a "$A416/system" "$CCUR/"
cp "$TEAM/dicts/4node/controlDict" "$CCUR/system/controlDict"     # team file, endTime 40
cp -a "$A416"/constant/{transportProperties,turbulenceProperties} "$CCUR/constant/"
ln -s "$MESH" "$CCUR/constant/polyMesh"
for p in "$A416"/processors416*; do ln -s "$p" "$CCUR/${p##*/}"; done
cp -a "$TIMER_SRC" "$CCUR/dynamicCode"

# --- w05: parallel prep (site module, like the team chain), reused if already complete ---
grep -q '^End' "$CW05/logFiles/10_decomposePar.log"
grep -q '^nDomainWeights (1 0.5);' "$CW05/system/include/caseDefinition"
prep_done=1
for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CW05/logFiles/$l.log" 2>/dev/null || prep_done=0; done
if [[ $prep_done == 1 ]]; then echo "REUSING prepared w05 case $CW05"; else
(
    set +eu
    module purge
    module load OpenFOAM/2512
    set -eu
    cd "$CW05"
    par=(mpirun -np 416 --map-by ppr:104:node --bind-to core)
    "${par[@]}" renumberMesh -constant -overwrite -parallel "${fhw05[@]}" > logFiles/20_renumberMesh.log 2>&1
    "${par[@]}" potentialFoam -initialiseUBCs -parallel "${fhw05[@]}" > logFiles/30_potentialFoam.log 2>&1
    "${par[@]}" applyBoundaryLayer -ybl "0.0450244" -parallel "${fhw05[@]}" > logFiles/40_applyBoundaryLayer.log 2>&1
    for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' logFiles/$l.log; done
)
fi
[[ -d $CW05/dynamicCode ]] || cp -a "$TIMER_SRC" "$CW05/dynamicCode"
echo "PREPARED_W05 $(grep -E 'Max number of processor patches' "$CW05/logFiles/10_decomposePar.log")"
echo "REFERENCE_CUR $(grep -E 'Max number of processor patches' "$A416/logFiles/10_decomposePar.log")"

# --- solver environment (team runtime settings) -------------------------------------
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
echo "simpleFoam=$(command -v simpleFoam)" > build-metadata.txt
mpirun -np 416 -ppn 104 hostname | sort | uniq -c > placement.txt
if [[ $(wc -l < placement.txt) != 4 ]] || awk '$1 != 104 {bad=1} END {exit !bad}' placement.txt; then
    echo "PLACEMENT_WRONG:"; cat placement.txt; exit 1
fi

run() {   # $1 arm (cur|w05)
    local arm=$1 cdir dir fh
    if [[ $arm == cur ]]; then cdir=$CCUR fh=("${fhcur[@]}"); else cdir=$CW05 fh=("${fhw05[@]}"); fi
    dir=$ROOT/run_$arm
    mkdir "$dir"
    cd "$cdir"
    [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == 40 ]]
    echo "RUN_START arm=$arm utc=$(date -u +%FT%TZ)"
    local status=0
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == 40 ]]
    grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS arm=$arm Cd_step40=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

run cur
run w05

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
cur, w05 = t('cur'), t('w05')
lines = ['# scotch domainWeights (1 0.5) vs current (1 0.01), 4 nodes x 104 ranks, 65M, 40 steps', '',
         'ICX OptSPR + Intel MPI 2021.17.2, scatter pinning, 4-node tuned Allreduce, glibc. One run each, same job.', '',
         '| Decomposition | Worst neighbours | Worst proc faces | s/step | vs current | Cd step 40 (0.43-0.44) |',
         '|---|---:|---:|---:|---:|---:|',
         f'| current (1 0.01) | 58 | 25,206 | {cur:.4f} | - | {cd("cur"):.4f} |',
         f'| (1 0.5) | 43 | 23,868 | {w05:.4f} | {100*(w05/cur-1):+.2f}% | {cd("w05"):.4f} |']
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "DW05_COMPLETE root=$ROOT"
