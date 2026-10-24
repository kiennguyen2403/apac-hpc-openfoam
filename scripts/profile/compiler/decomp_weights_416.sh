#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=40,mem=190gb,jobfs=10gb
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N decomp_wts_416
#PBS -j oe

# Score wall-weighted scotch decompositions (416 = 4 x 104, domains (4 104), link costs (2 1), as
# case416_w05) with decompScore before any decomposePar/solver run. Weights per cell
#   w = A*(internal faces) + B*(wall faces) + C, dithered integers (cellWeights -mode linear).
# Cost model from the imbalance diagnosis of the current setup (job 182137769, FIT.md): per-rank compute
# ~ 2.53 us/internal face + 5.23 us/proc face + 4.10 us/wall face + 1.29 ms/neighbour, i.e. a wall face
# ~3.2x an internal-face incidence. Candidates: unweighted (current), internal faces only, and wall
# weighting at 25/50/100/150% of the fitted ratio (A=10, B=8/16/32/48).
set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh
T=/scratch/ya41/kn3896/openfoam/tools/decompScore
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/decomp-weights/score_416_${PBS_JOBID%%.*}
W=$PBS_JOBFS/work
mkdir -p "$ROOT" "$W/case/constant" "$W/case/system"
exec > >(tee "$ROOT/job.log") 2>&1
ln -s "$MESH/constant/polyMesh" "$W/case/constant/polyMesh"
printf 'FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }\napplication none; startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1; writeControl timeStep; writeInterval 1;\n' > "$W/case/system/controlDict"
printf 'FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }\nnumberOfSubdomains 416;\nmethod scotch;\ncoeffs { domains (4 104); domainWeights (2 1); }\n' > "$W/w2"
cat > "$W/cands.txt" <<CANDS
unweighted $W/w2
faces $W/w2 10 0 0
wall25 $W/w2 10 8 0
wall50 $W/w2 10 16 0
wall100 $W/w2 10 32 0
wall150 $W/w2 10 48 0
CANDS
cp "$W/cands.txt" "$ROOT/"
set +eu
module purge; module load openmpi/4.1.7
source /apps/OpenFOAM/2512/etc/bashrc > /dev/null 2>&1
set -eu
"$T/bin/decompScore" -case "$W/case" -candidates "$W/cands.txt" -outDir "$ROOT" -nParallel 6 > "$ROOT/decompScore.log" 2>&1
grep -E '^CANDIDATE|^Candidates' "$ROOT/decompScore.log"
if cmp -s <(grep -v '^#' "$ROOT/unweighted.ranks" | cut -d' ' -f1-6) <(grep -v '^#' "$T/train_case416_w05.ranks"); then
    echo "CHECK_UNWEIGHTED_REPRODUCES_CASE416_W05: yes"
else
    echo "CHECK_UNWEIGHTED_REPRODUCES_CASE416_W05: NO"
fi
echo "DECOMP_WEIGHTS_COMPLETE root=$ROOT"
