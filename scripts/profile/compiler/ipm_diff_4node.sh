#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N ipm_diff_4n
#PBS -j oe

# IPM profile of the DEFAULT Gadi module (module load OpenFOAM/2512: GCC 8.5, Open MPI 4.1.7) on
# 4 nodes x 104 ranks, occDrivAerStaticMesh 65M, scotch domains (4 104) domainWeights (1 0.5).
# Reuses the prepared case416_w05 (decomposePar, restore0Dir, renumberMesh, potentialFoam,
# applyBoundaryLayer done with this same site module) -- no new decomposition.
# IPM as in the team's runs: module load ipm (Gadi's Open MPI picks it up), plain mpirun -np 416.
# Two runs in the same job, 5 and 25 steps: per-rank (wall - MPI) difference = cost of 20 steps with
# mesh read / startup cancelled. Training data for the per-step per-rank model (decomposition search;
# the whole-run fit from job 180805820 included ~47 s of startup and mispredicted w5 / sock211).

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
C=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/ipm-site/ipm_diff_${PBS_JOBID%%.*}
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "IPM_FAILED line=$LINENO status=$?"' ERR
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$C/logFiles/$l.log"; done
grep -q '^nDomainWeights (1 0.5);' "$C/system/include/caseDefinition"
cmp "$C/system/fvSolution" "$C/system/fvSolution.fixedIter"
ls "$C/dynamicCode/platforms/linux64GccDPInt32Opt/lib/"libwallClockTimeStatistics_*.so > /dev/null
if [[ -d $C/postProcessing ]]; then mv "$C/postProcessing" "$ROOT/postProcessing.before"; fi

set +eu
module purge
module load OpenFOAM/2512
module load ipm
set -eu
export WM_COMPILER_TYPE=system WM_CC=gcc WM_CXX=g++
{
    echo "simpleFoam=$(command -v simpleFoam)  WM_OPTIONS=$WM_OPTIONS"
    echo "mpirun=$(command -v mpirun)  $(mpirun --version 2>&1 | head -1)"
    echo "IPM_ROOT=${IPM_ROOT:-unset} IPM_LOG=${IPM_LOG:-unset}"
    module list 2>&1
} > "$ROOT/environment.txt"
export IPM_LOGDIR=$ROOT IPM_REPORT=full

run() {   # $1 steps
    local n=$1 d=$ROOT/steps$1 status=0
    mkdir -p "$d"
    export IPM_LOGDIR=$d
    cd "$C"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $n;/" system/controlDict
    [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == $n ]]
    echo "RUN_START steps=$n utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 900s mpirun -np 416 simpleFoam -parallel "${fh[@]}" > "$d/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$d/"; fi
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$d/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$d/simpleFoam.log") == $n ]]
    ls "$d"/*.ipm.xml > /dev/null
    echo "RUN_SUCCESS steps=$n $(grep 'Average wall-clock time per time step' "$d/simpleFoam.log" | tail -1)"
}
run 5
run 25
echo "IPM_COMPLETE root=$ROOT"
