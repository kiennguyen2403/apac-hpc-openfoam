#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N mesh_motorbike_M
#PBS -j oe

# HPC_motorbike Medium (~17M cells) for the one-node compiler comparison.
# Upstream recipe (committees/hpc develop 84c2624, Medium/v1912/AllmeshM) with
# the site OpenFOAM/2512 module, then a 104-way decomposition with initial fields.
# The earlier copy was decomposed after snappyHexMesh was OOM-killed: 20000 cells.

module purge
module load OpenFOAM/2512

src=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/HPC_motorbike/Medium/v1912
case=/scratch/ya41/kn3896/openfoam/contest-results/compiler-comparison/HPC_motorbike

[[ ! -e $case ]] || { echo "ERROR: $case already exists" >&2; exit 1; }
mkdir -p "$case"
cp -a "$src"/. "$case"/
cd "$case"
git -C "$src" rev-parse HEAD > source-commit.txt

. "$WM_PROJECT_DIR/bin/tools/RunFunctions"
set -e

# snappyHexMesh and the final decomposition both use all 104 cores.
foamDictionary -entry numberOfSubdomains -set 104 system/decomposeParDict

runApplication blockMesh
runApplication decomposePar
runParallel snappyHexMesh -overwrite
grep -q '^Finished meshing' log.snappyHexMesh || { echo "ERROR: snappyHexMesh did not finish" >&2; exit 1; }
runApplication reconstructParMesh -constant
runApplication reconstructPar -constant
foamDictionary -entry "pointAndNormalDict.basePoint" -set "( 0 4 0 )" system/mirrorMeshDict
runApplication mirrorMesh -overwrite
runApplication checkMesh

cells=$(awk '/^ *cells:/{print $2; exit}' log.checkMesh)
echo "checkMesh cells=$cells"
(( cells > 15000000 )) || { echo "ERROR: expected ~17M cells, got $cells" >&2; exit 1; }
grep -q 'Mesh OK' log.checkMesh || echo "WARNING: checkMesh did not report 'Mesh OK' (see log.checkMesh)"

# Replace the snappyHexMesh processor directories with the final decomposition plus 0/ fields.
cp -r 0.org 0
runApplication -s final decomposePar -force

nproc=$(find . -maxdepth 1 -type d -name 'processor[0-9]*' | wc -l)
[[ $nproc == 104 ]] || { echo "ERROR: $nproc processor directories" >&2; exit 1; }
total=0
for owner in processor*/constant/polyMesh/owner*; do
    n=$(zcat -f "$owner" | head -30 | grep -o 'nCells:[0-9]*' | cut -d: -f2)
    total=$((total + n))
done
echo "decomposed cells=$total"
(( total == cells )) || { echo "ERROR: decomposed $total cells, checkMesh $cells" >&2; exit 1; }
echo "MESH_SUCCESS cells=$cells processors=$nproc case=$case"
