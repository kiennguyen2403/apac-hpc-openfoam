#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=26,mem=128gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N numaflat_decomp
#PBS -j oe

# Topology-aware decomposition, step 1 (serial, small job): occDrivAerStaticMesh 65M into
# 416 = 4 nodes x 2 sockets x 4 NUMA x 13 cores, scotch domains (4 2 4 13), domainWeights
# (1 0.5 0.5 0.5) ("numaflat": best balance in the dry-run screen, job 180812660). Meant for
# COMPACT pinning (rank i on core i; NUMA k = cores 13k..13k+12 on normalsr). Collated output
# grouped per node (-ioRanks), then restore0Dir. Parallel prep runs in the 4-node job.

set -euo pipefail
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen
C=$ROOT/case416_numaflat
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

# An existing case408 is reused only if it holds no decomposition yet (a failed first attempt
# left one with copied dictionaries only); files below are overwritten.
if compgen -G "$C/processors416*" > /dev/null; then echo "ERROR: $C already decomposed" >&2; exit 1; fi
if [[ -L $C/constant/polyMesh ]]; then rm -- "$C/constant/polyMesh"; fi
mkdir -p "$C/constant" "$C/logFiles"
cp -a "$CASE_SRC/0.orig" "$CASE_SRC/system" "$C/"   # merges over an earlier copy
cp -a "$CASE_SRC"/constant/{transportProperties,turbulenceProperties} "$C/constant/"
ln -s "$MESH" "$C/constant/polyMesh"
cp "$TEAM/dicts/4node/caseDefinition"   "$C/system/include/caseDefinition"
cp "$TEAM/dicts/4node/decomposeParDict" "$C/system/decomposeParDict"
cp "$TEAM/dicts/4node/controlDict"      "$C/system/controlDict"   # team file, endTime 40 (screening)
sed -i -e 's/^nDomains .*/nDomains (4 2 4 13);/' -e 's/^nDomainWeights .*/nDomainWeights (1 0.5 0.5 0.5);/' "$C/system/include/caseDefinition"
grep -E '^(nCores|decompositionMethod|nDomains|nDomainWeights)' "$C/system/include/caseDefinition"

set +eu
module purge
module load OpenFOAM/2512
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
set -eu
cd "$C"
[[ $(foamDictionary -expand system/decomposeParDict | awk '/^numberOfSubdomains/ {sub(";","",$2); print $2}') == 416 ]]
decomposePar -constant "${fh[@]}" > logFiles/10_decomposePar.log 2>&1
grep -q '^End' logFiles/10_decomposePar.log
restore0Dir -processor > logFiles/15_restore0Dir.log 2>&1
ls -d processors416*
grep -E 'Max number of (cells|processor patches|faces between)' logFiles/10_decomposePar.log
echo "DECOMPOSE_NUMAFLAT_SUCCESS case=$C"
