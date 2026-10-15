#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N glibct_4n
#PBS -j oe

# glibc malloc tuning at 4 nodes: plain glibc vs glibcT (keep freed memory instead of returning it
# to the OS; no huge pages, so independent of node memory fragmentation). 1-node results: -0.4% and
# -0.7% (jobs 180771807, 180774706). 4 nodes x 104 ranks, occDrivAerStaticMesh 65M, 40 steps,
# prepared scotch (4 104)/(1 0.5) case (no prep), ICX OptSPR + Intel MPI, team scatter pinning,
# 4-node tuned Allreduce. One run each, same job/nodes.
#   glibc  : defaults
#   glibcT : MALLOC_MMAP_THRESHOLD_=32 MiB, MALLOC_TRIM_THRESHOLD_=4 GiB, MALLOC_TOP_PAD_=256 MiB

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
C=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/alloc-thp/glibct_4node_${PBS_JOBID%%.*}
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')
STEPS=40

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "GLIBCT_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$C/logFiles/$l.log"; done
grep -q '^nDomainWeights (1 0.5);' "$C/system/include/caseDefinition"
cmp "$C/system/fvSolution" "$C/system/fvSolution.fixedIter"
ls "$C/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
if [[ -d $C/postProcessing ]]; then mv "$C/postProcessing" "$ROOT/postProcessing.before"; fi

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

run() {   # $1 variant (glibc|glibcT)
    local v=$1 dir=$ROOT/run_$1 genv=() status=0
    [[ $v == glibcT ]] && genv=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)
    mkdir "$dir"
    cd "$C"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
    [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == $STEPS ]]
    echo "RUN_START variant=$v utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 "${genv[@]}" simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS variant=$v Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

run glibc
run glibcT
cp "$TEAM/dicts/4node/controlDict" "$C/system/controlDict"

python3 - "$ROOT" <<'EOF'
import re, sys
from pathlib import Path
root = Path(sys.argv[1])
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
def t(v): return float(pat.findall((root / f'run_{v}' / 'simpleFoam.log').read_text())[-1])
def cd(v):
    rows = [l.split() for p in (root / f'run_{v}' / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return float(rows[-1][1])
a, b = t('glibc'), t('glibcT')
lines = ['# glibc vs glibcT, 4 nodes x 104, occDrivAerStaticMesh 65M, scotch (1 0.5), 40 steps', '',
         'ICX OptSPR + Intel MPI, scatter pinning, 4-node tuned Allreduce. One run each, same job/nodes.', '',
         '| Allocator | s/step | vs glibc | Cd step 40 |', '|---|---:|---:|---:|',
         f'| glibc (default) | {a:.4f} | - | {cd("glibc"):.5f} |',
         f'| glibcT (MALLOC_MMAP_THRESHOLD_=32M, TRIM_THRESHOLD_=4G, TOP_PAD_=256M) | {b:.4f} | {100*(b/a-1):+.2f}% | {cd("glibcT"):.5f} |']
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "GLIBCT_COMPLETE root=$ROOT"
