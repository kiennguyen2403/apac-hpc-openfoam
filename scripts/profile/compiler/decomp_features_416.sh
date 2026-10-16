#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=12,mem=60gb,jobfs=1gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N decomp_feat
#PBS -j oe

# Per-rank features of the current decomposition (scotch (4 104)/(1 0.5) == link costs (2 1), i.e.
# case416_w05) incl. internal faces, wall faces and wall cells, for fitting a per-cell cost model
# against the per-rank times of the imbalance diagnosis (job 181350543). decompScore, site build.
set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh
T=/scratch/ya41/kn3896/openfoam/tools/decompScore
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/diag/features_416_${PBS_JOBID%%.*}
W=$PBS_JOBFS/work
mkdir -p "$ROOT" "$W/case/constant" "$W/case/system"
exec > >(tee "$ROOT/job.log") 2>&1
ln -s "$MESH/constant/polyMesh" "$W/case/constant/polyMesh"
printf 'FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }\napplication none; startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1; writeControl timeStep; writeInterval 1;\n' > "$W/case/system/controlDict"
printf 'FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }\nnumberOfSubdomains 416;\nmethod scotch;\ncoeffs { domains (4 104); domainWeights (2 1); }\n' > "$W/w2"
printf 'w2 %s/w2\n' "$W" > "$W/cands.txt"
set +eu
module purge; module load openmpi/4.1.7
source /apps/OpenFOAM/2512/etc/bashrc > /dev/null 2>&1
set -eu
"$T/bin/decompScore" -case "$W/case" -candidates "$W/cands.txt" -outDir "$ROOT" > "$ROOT/decompScore.log" 2>&1
grep -E '^CANDIDATE|^Candidates' "$ROOT/decompScore.log"
if cmp -s <(grep -v '^#' "$ROOT/w2.ranks" | cut -d' ' -f1-6) <(grep -v '^#' "$T/train_case416_w05.ranks"); then
    echo "CHECK_REPRODUCES_CASE416_W05: yes"
else
    echo "CHECK_REPRODUCES_CASE416_W05: NO"
fi
echo "FEATURES_COMPLETE root=$ROOT"
