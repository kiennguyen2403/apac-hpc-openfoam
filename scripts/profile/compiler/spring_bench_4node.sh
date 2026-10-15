#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N spring_4n
#PBS -j oe

# BONUS (experimental): renumberMesh twice -- default CuthillMcKee, then spring -- on the
# scotch domainWeights (1 0.5) 416-rank decomposition of occDrivAerStaticMesh 65M.
# Rules: changing renumberMesh method/dict is Bonus; Basic chain stays CuthillMcKee only.
# Dry-run screen (job 180787768): spring on top of CuthillMcKee cut the profile 5.44e11 -> 2.30e11.
#   cm     : case416_w05 (prepared: decomposePar, restore0Dir, renumberMesh CM, potentialFoam,
#            applyBoundaryLayer) -- reference, untouched
#   spring : copy of the same case + renumberMesh (spring) + restore0Dir -processor +
#            potentialFoam -initialiseUBCs + applyBoundaryLayer -ybl 0.0450244
# One 40-step run each, same job/nodes, ICX OptSPR + Intel MPI, team scatter pinning, 4-node
# tuned Allreduce, glibc. Bonus check F3: Cd within +-0.5% of the Basic baseline.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
CCM=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
CSP=/scratch/ya41/kn3896/openfoam/contest-results/renumber-screen/case416_w05_spring
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/renumber-screen/bench_spring_${PBS_JOBID%%.*}
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "SPRING_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for c in "$CCM" "$CSP"; do
    for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$c/logFiles/$l.log"; done
    grep -q '^nDomainWeights (1 0.5);' "$c/system/include/caseDefinition"
    ls "$c/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
done
printf 'FoamFile { version 2.0; format ascii; class dictionary; object renumberMeshDict; }\nmethod spring;\nspringCoeffs { maxCo 0.01; freezeFraction 0.999; maxIter 1000; }\n' > "$ROOT/renumberMeshDict.spring"

# --- spring: second renumber pass + rebuild fields as in the chain (site module) ----------
if grep -q '^End' "$CSP/logFiles/41_applyBoundaryLayer_spring.log" 2>/dev/null; then
    echo "REUSING spring-prepared case $CSP"
else
(
    set +eu
    module purge
    module load OpenFOAM/2512
    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"
    set -eu
    cd "$CSP"
    par=(mpirun -np 416 --map-by ppr:104:node --bind-to core)
    "${par[@]}" renumberMesh -constant -overwrite -parallel "${fh[@]}" -dict "$ROOT/renumberMeshDict.spring" \
        > logFiles/21_renumberMesh_spring.log 2>&1
    grep -q '^End' logFiles/21_renumberMesh_spring.log
    grep -q 'renumber-method: spring' logFiles/21_renumberMesh_spring.log
    restore0Dir -processor > logFiles/22_restore0Dir_spring.log 2>&1
    "${par[@]}" potentialFoam -initialiseUBCs -parallel "${fh[@]}" > logFiles/31_potentialFoam_spring.log 2>&1
    "${par[@]}" applyBoundaryLayer -ybl "0.0450244" -parallel "${fh[@]}" > logFiles/41_applyBoundaryLayer_spring.log 2>&1
    for l in 31_potentialFoam_spring 41_applyBoundaryLayer_spring; do grep -q '^End' logFiles/$l.log; done
)
fi
echo "PREPARED_SPRING $(awk '/^After renumbering/{f=1} f && /band|profile/ {printf "%s=%s ", $1, $NF}' "$CSP/logFiles/21_renumberMesh_spring.log")"
echo "REFERENCE_CM    $(awk '/^After renumbering/{f=1} f && /band|profile/ {printf "%s=%s ", $1, $NF}' "$CCM/logFiles/20_renumberMesh.log")"

# --- solver environment (team runtime settings) -----------------------------------------
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
mpirun -np 416 -ppn 104 hostname | sort | uniq -c > placement.txt
if [[ $(wc -l < placement.txt) != 4 ]] || awk '$1 != 104 {bad=1} END {exit !bad}' placement.txt; then
    echo "PLACEMENT_WRONG:"; cat placement.txt; exit 1
fi

run() {   # $1 arm (cm|spring)
    local arm=$1 cdir dir
    if [[ $arm == cm ]]; then cdir=$CCM; else cdir=$CSP; fi
    dir=$ROOT/run_$arm
    mkdir "$dir"
    cd "$cdir"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
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

run cm
run spring

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
cm, sp = t('cm'), t('spring')
dcd = 100 * (cd('spring') / cd('cm') - 1)
lines = ['# BONUS: renumberMesh CuthillMcKee + spring vs CuthillMcKee only, 4 nodes x 104, 65M, 40 steps', '',
         'scotch domainWeights (1 0.5); ICX OptSPR + Intel MPI, scatter pinning, tuned Allreduce, glibc. One run each, same job.', '',
         '| Ordering | s/step | vs CM | Cd step 40 | Cd vs CM (F3: within +-0.5%) |', '|---|---:|---:|---:|---:|',
         f'| CuthillMcKee (Basic) | {cm:.4f} | - | {cd("cm"):.5f} | - |',
         f'| CuthillMcKee + spring (Bonus) | {sp:.4f} | {100*(sp/cm-1):+.2f}% | {cd("spring"):.5f} | {dcd:+.3f}% {"PASS" if abs(dcd) <= 0.5 else "FAIL"} |']
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "SPRING_COMPLETE root=$ROOT"
