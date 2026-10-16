#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=40,mem=190gb,jobfs=2gb
#PBS -l walltime=00:50:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N weighted_dec
#PBS -j oe

# BONUS (weightField): weighted scotch decompositions for the 4-node case, to balance estimated
# work instead of cell count (imbalance diagnosis job 181350543: per-rank cost follows the
# subdomain; fit job 181351948: cost ~ faces, +wall faces). Team case + team 4-node dicts
# (scotch (4 104)/(1 0.5)) unchanged except "weightField cellWeight;" in decomposeParDict.
#   wface : cellWeights -mode faces     (w = 4*internal faces + 7*wall faces per cell)
#   wmeas : cellWeights -mode measured  (faces weight scaled per rank by measured cost)
# Then decomposePar (collated per node) and restore0Dir; parallel prep runs in the 4-node job.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
CW=/scratch/ya41/kn3896/openfoam/tools/cellWeights/bin/cellWeights
W=/scratch/ya41/kn3896/openfoam/contest-results/weighted
COSTS=$W/costs_181350543.txt
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')
exec > >(tee "$W/decompose_${PBS_JOBID%%.*}.log") 2>&1
trap 'echo "WEIGHTED_DECOMPOSE_FAILED line=$LINENO status=$?"' ERR

# decomposition that the measured costs belong to (== case416_w05, verified job 181351948)
printf 'FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }\nnumberOfSubdomains 416;\nmethod scotch;\ncoeffs { domains (4 104); domainWeights (2 1); }\n' > "$W/dict_w2_measured"

setup() {   # $1 name
    local C=$W/case416_$1
    if [[ -e $C ]]; then echo "ERROR: $C exists" >&2; exit 1; fi
    mkdir -p "$C/constant" "$C/logFiles" "$C/0"
    cp -a "$CASE_SRC/0.orig" "$CASE_SRC/system" "$C/"
    cp -a "$CASE_SRC"/constant/{transportProperties,turbulenceProperties} "$C/constant/"
    ln -s "$MESH" "$C/constant/polyMesh"
    cp "$TEAM/dicts/4node/caseDefinition"   "$C/system/include/caseDefinition"
    cp "$TEAM/dicts/4node/decomposeParDict" "$C/system/decomposeParDict"
    cp "$TEAM/dicts/4node/controlDict"      "$C/system/controlDict"
    printf '\nweightField cellWeight;\n' >> "$C/system/decomposeParDict"
    grep -q '^weightField cellWeight;' "$C/system/decomposeParDict"
}
setup wface
setup wmeas

set +eu
module purge
module load OpenFOAM/2512
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
set -eu
type restore0Dir > /dev/null

echo "WEIGHTS_START utc=$(date -u +%FT%TZ)"
"$CW" -case "$W/case416_wface" -mode faces > "$W/case416_wface/logFiles/05_cellWeights.log" 2>&1 &
p1=$!
"$CW" -case "$W/case416_wmeas" -mode measured -dict "$W/dict_w2_measured" -costs "$COSTS" > "$W/case416_wmeas/logFiles/05_cellWeights.log" 2>&1 &
p2=$!
wait $p1; wait $p2
for n in wface wmeas; do grep -E 'measured:|weights:' "$W/case416_$n/logFiles/05_cellWeights.log" | sed "s/^/$n: /"; ls -la "$W/case416_$n/0/cellWeight"; done

echo "DECOMPOSE_START utc=$(date -u +%FT%TZ)"
pids=()
for n in wface wmeas; do
    (
        cd "$W/case416_$n"
        [[ $(foamDictionary -expand system/decomposeParDict | awk '/^numberOfSubdomains/ {sub(";","",$2); print $2}') == 416 ]]
        decomposePar -constant "${fh[@]}" > logFiles/10_decomposePar.log 2>&1
        grep -q '^End' logFiles/10_decomposePar.log
        restore0Dir -processor > logFiles/15_restore0Dir.log 2>&1
    ) &
    pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
for n in wface wmeas; do
    echo "$n: $(grep -E 'Max number of (cells|processor patches|faces between)' "$W/case416_$n/logFiles/10_decomposePar.log" | tr -s ' ' | tr '\n' ' ')"
    ls "$W/case416_$n/processors416_0-103/0" | tr '\n' ' '; echo
done
echo "WEIGHTED_DECOMPOSE_SUCCESS"
