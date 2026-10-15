#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=26,mem=128gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N final_dec208
#PBS -j oe

# Final runs, 2-node case: decompose occDrivAerStaticMesh 65M into 208 = 2 x 104 parts with the
# team's 2-node dicts unchanged (scotch, nDomains (2 104), nDomainWeights (1 0.5)), collated per
# node (-ioRanks (0 104)), then restore0Dir. The parallel prep runs in the 2-node final job.

set -euo pipefail
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/final
C=$ROOT/case208
fh=(-fileHandler collated -ioRanks '(0 104)')

if [[ -e $C ]]; then echo "ERROR: $C exists" >&2; exit 1; fi
mkdir -p "$C/constant" "$C/logFiles"
cp -a "$CASE_SRC/0.orig" "$CASE_SRC/system" "$C/"   # merges over an earlier copy
cp -a "$CASE_SRC"/constant/{transportProperties,turbulenceProperties} "$C/constant/"
ln -s "$MESH" "$C/constant/polyMesh"
cp "$TEAM/dicts/2node/caseDefinition"   "$C/system/include/caseDefinition"
cp "$TEAM/dicts/2node/decomposeParDict" "$C/system/decomposeParDict"
cp "$TEAM/dicts/2node/controlDict"      "$C/system/controlDict"   # team file, endTime 80
grep -E '^(nCores|decompositionMethod|nDomains|nDomainWeights)' "$C/system/include/caseDefinition"

set +eu
module purge
module load OpenFOAM/2512
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
set -eu
cd "$C"
[[ $(foamDictionary -expand system/decomposeParDict | awk '/^numberOfSubdomains/ {sub(";","",$2); print $2}') == 208 ]]
decomposePar -constant "${fh[@]}" > logFiles/10_decomposePar.log 2>&1
grep -q '^End' logFiles/10_decomposePar.log
restore0Dir -processor > logFiles/15_restore0Dir.log 2>&1
ls -d processors208*
grep -E 'Max number of (cells|processor patches|faces between)' logFiles/10_decomposePar.log
echo "DECOMPOSE_208_SUCCESS case=$C"
