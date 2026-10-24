#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=40,mem=190gb,jobfs=10gb
#PBS -l walltime=03:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N kahip_score_416
#PBS -j oe

# Score KaHIP decompositions for 416 ranks (4 nodes x 104) with decompScore (per-rank cells, processor
# faces, neighbours, wall faces; no processor dirs written), to tune the team's KaHIP partition
# (multiLevel: nodes kahip strong, cores kahip eco, default imbalance 1%: largest rank +2.0% cells, the
# compounded tolerance of two levels; job 182167199). References:
#   scotch_ref - scotch (4 104)/(2 1) = case416_w05 (must reproduce tools/decompScore/train_case416_w05.ranks)
#   team       - the team's multiLevel KaHIP (must reproduce the team partition's statistics)
# Variants: imbalance 0.3% / 0.1% per level, strong/strong, fast/fast, flat 416-way eco/fast, and cell
# weights w = A*(internal faces) + B*(wall faces) + C (as cellWeights -mode linear) on strong/eco 0.3%.
set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh
T=/scratch/ya41/kn3896/openfoam/tools/decompScore
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/kahip/score_416_${PBS_JOBID%%.*}
W=$PBS_JOBFS/work
mkdir -p "$ROOT" "$W/case/constant" "$W/case/system" "$W/dicts"
exec > >(tee "$ROOT/job.log") 2>&1
ln -s "$MESH/constant/polyMesh" "$W/case/constant/polyMesh"
printf 'FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }\napplication none; startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1; writeControl timeStep; writeInterval 1;\n' > "$W/case/system/controlDict"

hdr='FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }'
printf '%s\nnumberOfSubdomains 416;\nmethod scotch;\ncoeffs { domains (4 104); domainWeights (2 1); }\n' "$hdr" > "$W/dicts/scotch_ref"
ml() {   # $1 name, $2 nodes config, $3 cores config, $4 imbalance ("" = KaHIP default)
    local imb=""; [[ -n $4 ]] && imb="imbalance $4;"
    printf '%s\nnumberOfSubdomains 416;\nmethod multiLevel;\nmultiLevelCoeffs\n{\n    nodes { numberOfSubdomains 4; method kahip; coeffs { config %s; %s } }\n    cores { numberOfSubdomains 104; method kahip; coeffs { config %s; %s } }\n}\n' \
        "$hdr" "$2" "$imb" "$3" "$imb" > "$W/dicts/$1"
}
flat() {   # $1 name, $2 config, $3 imbalance
    printf '%s\nnumberOfSubdomains 416;\nmethod kahip;\ncoeffs { config %s; imbalance %s; }\n' "$hdr" "$2" "$3" > "$W/dicts/$1"
}
ml team strong eco ""
ml ml_se_i003 strong eco 0.003
ml ml_se_i001 strong eco 0.001
ml ml_ss_i003 strong strong 0.003
ml ml_ff_i003 fast fast 0.003
flat flat_eco_i003 eco 0.003
flat flat_fast_i003 fast 0.003
cat > "$W/cands.txt" <<CANDS
scotch_ref $W/dicts/scotch_ref
team $W/dicts/team
ml_se_i003 $W/dicts/ml_se_i003
ml_se_i001 $W/dicts/ml_se_i001
ml_ss_i003 $W/dicts/ml_ss_i003
ml_ff_i003 $W/dicts/ml_ff_i003
flat_eco_i003 $W/dicts/flat_eco_i003
flat_fast_i003 $W/dicts/flat_fast_i003
ml_se_i003_wface $W/dicts/ml_se_i003 10 0 0
ml_se_i003_wwall $W/dicts/ml_se_i003 10 16 0
CANDS
cp -a "$W/dicts" "$W/cands.txt" "$ROOT/"
set +eu
module purge; module load openmpi/4.1.7
source /apps/OpenFOAM/2512/etc/bashrc > /dev/null 2>&1
set -eu
"$T/bin/decompScore" -case "$W/case" -candidates "$W/cands.txt" -outDir "$ROOT" -nParallel 6 > "$ROOT/decompScore.log" 2>&1
grep -E '^CANDIDATE|^Candidates' "$ROOT/decompScore.log"
if cmp -s <(grep -v '^#' "$ROOT/scotch_ref.ranks" | cut -d' ' -f1-6) <(grep -v '^#' "$T/train_case416_w05.ranks"); then
    echo "CHECK_SCOTCH_REF_REPRODUCES_CASE416_W05: yes"
else
    echo "CHECK_SCOTCH_REF_REPRODUCES_CASE416_W05: NO"
fi
echo "KAHIP_SCORE_COMPLETE root=$ROOT"
