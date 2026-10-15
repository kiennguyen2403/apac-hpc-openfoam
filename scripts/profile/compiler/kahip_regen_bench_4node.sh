#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N kahip_regen_4n
#PBS -j oe

# BONUS (KaHIP is not Basic-legal): the team's KaHIP partition (an0642 kahip_4_node: multiLevel,
# nodes kahip strong / cores kahip eco; max 37 neighbours, inter-node faces 140k vs 275k) vs the
# current scotch (4 104)/(1 0.5), with the current Bonus setup (fixedIterSolvers plugin + fused
# vector solve). The team case's processor meshes are intact (decomposed, renumbered on Gadi with
# v2512) but its fields were deleted, so:
#   - our case links each processors416_* constant/ to the team's (read-only, nothing written)
#   - fields are regenerated in our own 0/ with the team chain: restore0Dir -processor,
#     potentialFoam -initialiseUBCs, applyBoundaryLayer (site OpenFOAM/2512, hostCollated)
# Then paired, 4 nodes x 104, 10 steps each, same job. The scotch arm runs in a linked copy of
# case416_w05 (read-only processor directories).

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
KTEAM=/scratch/ya41/an0642/case/kahip_4_node
CUR=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/kahip/regen_bench_${PBS_JOBID%%.*}
CK=$ROOT/case_kahip
CS=$ROOT/case_scotch
STEPS=10
fh_scotch=(-fileHandler collated -ioRanks '(0 104 208 312)')
fh_kahip=(-fileHandler hostCollated)

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "KAHIP_REGEN_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh; do grep -q '^End' "$KTEAM/logFiles/$l.log"; done
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CUR/logFiles/$l.log"; done
[[ -f $LIB && -f $FLIB ]]
cp "$CUR/system/fvSolution.fixedIter" fvSolution.basic
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.basic > fvSolution.bonus
[[ $(grep -c 'smoothSolverFixed;' fvSolution.bonus) == 3 && $(grep -c 'PCGFixed;' fvSolution.bonus) == 1 ]]

setup_case() {   # $1 case dir: our case setup (system/, constant/ files, timer)
    mkdir -p "$1/constant" "$1/logFiles"
    cp -a "$CUR/system" "$1/"
    cp fvSolution.basic "$1/system/fvSolution"
    cp "$TEAM/dicts/4node/controlDict" "$1/system/controlDict"
    cp -a "$CUR"/constant/{transportProperties,turbulenceProperties} "$1/constant/"
    ln -s "$MESH" "$1/constant/polyMesh"
    cp -a "$CUR/dynamicCode" "$1/dynamicCode"
}

# --- scotch arm: links to the prepared case416_w05 processor directories --------------------
setup_case "$CS"
for p in "$CUR"/processors416_*; do ln -s "$p" "$CS/${p##*/}"; done

# --- KaHIP arm: team processor meshes (read-only links) + our own 0/ -------------------------
setup_case "$CK"
cp -a "$CASE_SRC/0.orig" "$CK/"
for p in "$KTEAM"/processors416_*; do
    d=$CK/${p##*/}
    mkdir "$d"
    ln -s "$p/constant" "$d/constant"
done
echo "KaHIP case: $(ls -d "$CK"/processors416_* | wc -l) processor groups, constant/ linked to $KTEAM"
grep -E 'Max number of (cells|processor patches|faces between)' "$KTEAM/logFiles/10_decomposePar.log"
(
    set +eu
    module purge
    module load OpenFOAM/2512
    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"
    set -eu
    cd "$CK"
    par=(mpirun -np 416 --map-by ppr:104:node --bind-to core)
    echo "PREP_START utc=$(date -u +%FT%TZ)"
    restore0Dir -processor > logFiles/15_restore0Dir.log 2>&1
    "${par[@]}" potentialFoam -initialiseUBCs -parallel "${fh_kahip[@]}" > logFiles/30_potentialFoam.log 2>&1
    "${par[@]}" applyBoundaryLayer -ybl "0.0450244" -parallel "${fh_kahip[@]}" > logFiles/40_applyBoundaryLayer.log 2>&1
    for l in 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' logFiles/$l.log; done
    echo "PREP_DONE utc=$(date -u +%FT%TZ) fields: $(ls processors416_0-103/0 | tr '\n' ' ')"
)

# --- solver environment (current Bonus setup) ------------------------------------------
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

run() {   # $1 arm (scotch|kahip)
    local a=$1 dir=$ROOT/run_$1 status=0 casedir
    if [[ $a == scotch ]]; then casedir=$CS; else casedir=$CK; fi
    local -n fh=fh_$a
    mkdir "$dir"
    cd "$casedir"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    echo "libs (\"$LIB\");" >> system/controlDict
    cp "$ROOT/fvSolution.bonus" system/fvSolution
    echo "RUN_START arm=$a utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 "${glibcT[@]}" -genv LD_PRELOAD "$FLIB" -genv FUSED_VECTOR_SOLVE_REPORT 1 \
        simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    grep -q 'fusedVectorSolve: U fused path' "$dir/simpleFoam.log"
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
run scotch
run kahip

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
s, k = t('scotch'), t('kahip')
out = [f'# Team KaHIP partition (fields regenerated) vs scotch (1 0.5), 4 nodes x 104, {steps} steps, Bonus setup, same job', '',
       '| Partition | max nbrs | inter-node faces | s/step | vs scotch | Cd step ' + steps + ' | Cd vs scotch |', '|---|---:|---:|---:|---:|---:|---:|',
       f'| scotch (4 104)/(1 0.5) | 43 | 275k | {s:.4f} | - | {cd("scotch"):.8f} | - |',
       f'| KaHIP multiLevel (team) | 37 | 140k | {k:.4f} | {100 * (k / s - 1):+.2f}% | {cd("kahip"):.8f} | {100 * (cd("kahip") / cd("scotch") - 1):+.4f}% |']
(root / 'REPORT.md').write_text('\n'.join(out) + '\n'); print('\n'.join(out))
PY
echo "KAHIP_REGEN_COMPLETE root=$ROOT"
