#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N aggl_4n
#PBS -j oe

# BONUS (experimental, F2): GAMG processorAgglomerator for p -- "fvSolution settings that move
# work between MPI ranks without changing it". 4 nodes x 104, occDrivAerStaticMesh 65M, 40 steps.
# fvSolution.fixedIter p: GAMG/GaussSeidel maxIter 12, coarsest level PCG/DIC 35 iterations across
# all 416 ranks (~10 cells/rank there: latency-bound). Variants (only this entry added to p):
#   base   : as shipped (no processorAgglomerator)
#   mc1    : processorAgglomerator masterCoarsest;              (coarsest level on one master)
#   mc104  : processorAgglomerator masterCoarsest; nProcessorsPerMaster 104;  (one master per node)
# Reuses the prepared scotch domainWeights (1 0.5) case (no prep). ICX OptSPR + Intel MPI, team
# scatter pinning, 4-node tuned Allreduce, glibc. Bonus check F3: Cd within +-0.5% of base.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
C=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/aggl/bench_${PBS_JOBID%%.*}
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')
VARIANTS="${AGGL_VARIANTS:-base mc1 mc104}"

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "AGGL_FAILED line=$LINENO status=$?"; cp "$C/system/fvSolution.fixedIter" "$C/system/fvSolution"' ERR
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

set_variant() {   # writes system/fvSolution for the variant (only p gets new entries)
    cp "$C/system/fvSolution.fixedIter" "$C/system/fvSolution"
    case $1 in
        base)  ;;
        mc1)   foamDictionary -entry solvers.p.processorAgglomerator -add masterCoarsest "$C/system/fvSolution" > /dev/null ;;
        mc104) foamDictionary -entry solvers.p.processorAgglomerator -add masterCoarsest "$C/system/fvSolution" > /dev/null
               foamDictionary -entry solvers.p.nProcessorsPerMaster -add 104 "$C/system/fvSolution" > /dev/null ;;
    esac
    foamDictionary -entry solvers.p "$C/system/fvSolution" > "$ROOT/run_$1/fvSolution.p.txt"
}

run() {   # $1 variant; a failed variant is recorded, the job continues
    local v=$1 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    set_variant "$v"
    cd "$C"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == 40 ]]
    echo "RUN_START variant=$v utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np 416 -ppn 104 simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cp "$C/system/fvSolution.fixedIter" "$C/system/fvSolution"
    cd "$ROOT"
    if [[ $status != 0 ]] || ! grep -qx End "$dir/simpleFoam.log" || ! grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log" \
       || [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") != 40 ]]; then
        echo "RUN_FAILED variant=$v status=$status: $(grep -m1 -iE 'FATAL|error|signal' "$dir/simpleFoam.log" | cut -c1-160)"
        return 0
    fi
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS variant=$v Cd_step40=$cd $(grep -m1 'GAMG:  Solving for p' "$dir/simpleFoam.log" | grep -oE 'Final residual = [0-9.e+-]+') $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

for v in $VARIANTS; do run "$v"; done
cmp "$C/system/fvSolution" "$C/system/fvSolution.fixedIter"   # case left as shipped

python3 - "$ROOT" $VARIANTS <<'EOF'
import re, sys
from pathlib import Path
root, variants = Path(sys.argv[1]), sys.argv[2:]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
def t(v):
    f = root / f'run_{v}' / 'simpleFoam.log'
    m = pat.findall(f.read_text()) if f.exists() else []
    return float(m[-1]) if m else None
def cd(v):
    rows = [l.split() for p in (root / f'run_{v}' / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return float(rows[-1][1]) if rows else None
base_t, base_cd = t('base'), cd('base')
lines = ['# BONUS: GAMG processorAgglomerator for p, 4 nodes x 104, 65M, 40 steps', '',
         'scotch domainWeights (1 0.5); ICX OptSPR + Intel MPI, scatter pinning, tuned Allreduce, glibc. One run each, same job.', '',
         '| Variant | s/step | vs base | Cd step 40 | Cd vs base (F3 +-0.5%) |', '|---|---:|---:|---:|---:|']
for v in variants:
    tv, cv = t(v), cd(v)
    if tv is None: lines.append(f'| {v} | failed | | | |'); continue
    rel = f'{100*(tv/base_t-1):+.2f}%' if base_t and v != 'base' else '-'
    dcd = f'{100*(cv/base_cd-1):+.3f}% {"PASS" if abs(cv/base_cd-1) <= 0.005 else "FAIL"}' if base_cd and v != 'base' else '-'
    lines.append(f'| {v} | {tv:.4f} | {rel} | {cv:.5f} | {dcd} |')
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "AGGL_COMPLETE root=$ROOT"
