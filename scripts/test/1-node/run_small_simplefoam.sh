#!/bin/bash -l
#PBS -N small_simplefoam
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=8GB
#PBS -l walltime=00:10:00
#PBS -l wd
#PBS -j oe
#PBS -o /home/565/kn3896/scratch/openfoam/logs/tests/small_simplefoam.out

set -eEo pipefail

ROOT=/home/565/kn3896/scratch/openfoam
FOAM_DIR="$ROOT/OpenFOAM-dev"
SOURCE_CASE="$FOAM_DIR/tutorials/incompressible/simpleFoam/pitzDaily"
CASE_DIR="$ROOT/tests/simpleFoam_pitzDaily_${PBS_JOBID:-manual}"
NPROCS="${PBS_NCPUS:-4}"

module purge
module load intel-compiler-llvm/2025.3.2 openmpi/5.0.8 metis/5.1.0 boost/1.84.0

export WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx
export WM_MPLIB=SYSTEMOPENMPI OMPI_CC=icx OMPI_CXX=icpx
source "$FOAM_DIR/etc/bashrc" || true
export WM_MPLIB=SYSTEMOPENMPI FOAM_MPI=sys-openmpi MPI_ARCH_PATH="$OPENMPI_ROOT"
export LD_LIBRARY_PATH="$FOAM_LIBBIN/$FOAM_MPI:$FOAM_LIBBIN:$OPENMPI_ROOT/lib:${LD_LIBRARY_PATH:-}"

mkdir -p "$ROOT/tests" "$ROOT/logs/tests"
cp -a "$SOURCE_CASE" "$CASE_DIR"

foamDictionary -entry endTime -set 100 "$CASE_DIR/system/controlDict"
foamDictionary -entry functions -remove "$CASE_DIR/system/controlDict"

cat > "$CASE_DIR/system/decomposeParDict" <<EOF
FoamFile
{
    format ascii;
    class dictionary;
    object decomposeParDict;
}
numberOfSubdomains $NPROCS;
method hierarchical;
hierarchicalCoeffs
{
    n (2 2 1);
    order xyz;
}
EOF

blockMesh -case "$CASE_DIR" > "$CASE_DIR/log.blockMesh" 2>&1
checkMesh -case "$CASE_DIR" > "$CASE_DIR/log.checkMesh" 2>&1
decomposePar -case "$CASE_DIR" > "$CASE_DIR/log.decomposePar" 2>&1

/usr/bin/time -f 'elapsed_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kb=%M' \
    -o "$CASE_DIR/time.simpleFoam.txt" \
    mpiexec -n "$NPROCS" simpleFoam -parallel -case "$CASE_DIR" \
    > "$CASE_DIR/log.simpleFoam" 2>&1

grep -q '^End$' "$CASE_DIR/log.simpleFoam"
! grep -Eq 'FOAM FATAL|Segmentation fault' "$CASE_DIR/log.simpleFoam"

echo "simpleFoam test completed successfully"
echo "Case: $CASE_DIR"
cat "$CASE_DIR/time.simpleFoam.txt"
grep -E '^(Time =|ExecutionTime =)' "$CASE_DIR/log.simpleFoam" | tail -4
