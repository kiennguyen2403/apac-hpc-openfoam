#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=26,mem=128gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N renum_raw_decomp
#PBS -j oe

# Option A, step 1: a raw (not renumbered) 104-rank decomposition of occDrivAerStaticMesh 65M with
# the team 1-node settings (plain scotch), so renumberMesh methods can be screened exactly as
# the mandated chain applies them (single pass on the decomposed mesh). Mesh only: no fields.

set -euo pipefail
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
C=/scratch/ya41/kn3896/openfoam/contest-results/renumber-screen/case104_raw

if compgen -G "$C/processors104*" > /dev/null; then echo "ERROR: $C already decomposed" >&2; exit 1; fi
mkdir -p "$C/constant" "$C/logFiles"
cp -a "$CASE_SRC/system" "$C/"
cp -a "$CASE_SRC"/constant/{transportProperties,turbulenceProperties} "$C/constant/"
ln -sfn "$MESH" "$C/constant/polyMesh"
cp "$TEAM/dicts/1node/caseDefinition"   "$C/system/include/caseDefinition"
cp "$TEAM/dicts/1node/decomposeParDict" "$C/system/decomposeParDict"
cp "$TEAM/dicts/1node/controlDict"      "$C/system/controlDict"

set +eu
module purge
module load OpenFOAM/2512
set -eu
cd "$C"
[[ $(foamDictionary -expand system/decomposeParDict | awk '/^numberOfSubdomains/ {sub(";","",$2); print $2}') == 104 ]]
decomposePar -constant -fileHandler collated > logFiles/10_decomposePar.log 2>&1
grep -q '^End' logFiles/10_decomposePar.log
ls -d processors104*
echo "DECOMPOSE_RAW104_SUCCESS case=$C"
