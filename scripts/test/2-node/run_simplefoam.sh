#!/bin/bash -l

#PBS -N simplefoam_2nodes
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=96
#PBS -l mem=192GB
#PBS -l walltime=00:15:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -j oe
#PBS -o /home/565/kn3896/scratch/openfoam/logs/tests/simplefoam_2nodes.out

set -eEo pipefail

ROOT=/home/565/kn3896/scratch/openfoam
FOAM_DIR="$ROOT/OpenFOAM-dev"
SOURCE_CASE="$FOAM_DIR/tutorials/incompressible/simpleFoam/pitzDaily"
CASE_DIR="$ROOT/tests/simpleFoam_pitzDaily_2nodes_${PBS_JOBID:-manual}"
NPROCS="${PBS_NCPUS:-96}"
mkdir -p "$ROOT/logs/tests"
DEBUG_LOG="$ROOT/logs/tests/simplefoam_2nodes_${PBS_JOBID:-manual}.log"
exec > >(tee "$DEBUG_LOG") 2>&1
trap 'rc=$?; echo "ERROR: line $LINENO: $BASH_COMMAND (exit $rc)" >&2' ERR

module purge
module load intel-compiler-llvm/2025.3.2
module load openmpi/5.0.8
module load metis/5.1.0
module load boost/1.84.0

export WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx
export WM_MPLIB=SYSTEMOPENMPI OMPI_CC=icx OMPI_CXX=icpx
source "$FOAM_DIR/etc/bashrc" \
    WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
    WM_MPLIB=SYSTEMOPENMPI || true

export WM_MPLIB=SYSTEMOPENMPI FOAM_MPI=sys-openmpi MPI_ARCH_PATH="$OPENMPI_ROOT"
export LD_LIBRARY_PATH="$FOAM_LIBBIN/$FOAM_MPI:$FOAM_LIBBIN:$OPENMPI_ROOT/lib:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=1

[[ "$NPROCS" == 96 ]] || {
    echo "ERROR: expected 96 allocated CPUs, got $NPROCS" >&2
    exit 1
}
[[ -r "${PBS_NODEFILE:-}" ]] || {
    echo "ERROR: PBS_NODEFILE is unavailable" >&2
    exit 1
}

NODE_COUNT="$(sort -u "$PBS_NODEFILE" | wc -l)"
[[ "$NODE_COUNT" == 2 ]] || {
    echo "ERROR: expected two nodes, got $NODE_COUNT" >&2
    sort "$PBS_NODEFILE" | uniq -c >&2
    exit 1
}

mkdir -p "$ROOT/tests" "$ROOT/logs/tests"
cp -a "$SOURCE_CASE" "$CASE_DIR"

# Short functional test only. pitzDaily is too small for performance scaling.
foamDictionary -entry endTime -set 25 "$CASE_DIR/system/controlDict"
foamDictionary -entry writeInterval -set 25 "$CASE_DIR/system/controlDict"
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
    n (12 8 1);
    order xyz;
}
EOF

echo "Allocated nodes:"
sort "$PBS_NODEFILE" | uniq -c
echo "MPI ranks: $NPROCS"

blockMesh -case "$CASE_DIR" > "$CASE_DIR/log.blockMesh" 2>&1
checkMesh -case "$CASE_DIR" > "$CASE_DIR/log.checkMesh" 2>&1
decomposePar -case "$CASE_DIR" > "$CASE_DIR/log.decomposePar" 2>&1

/usr/bin/time \
    -f 'elapsed_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kb=%M' \
    -o "$CASE_DIR/time.simpleFoam.txt" \
    mpiexec -n "$NPROCS" --map-by ppr:48:node:PE=1 --bind-to core \
        simpleFoam -parallel -case "$CASE_DIR" \
    > "$CASE_DIR/log.simpleFoam" 2>&1

grep -q '^End$' "$CASE_DIR/log.simpleFoam"
! grep -Eq 'FOAM FATAL|Segmentation fault' "$CASE_DIR/log.simpleFoam"

echo "Two-node simpleFoam test completed successfully"
echo "Case: $CASE_DIR"
cat "$CASE_DIR/time.simpleFoam.txt"
grep -E '^(Time =|ExecutionTime =)' "$CASE_DIR/log.simpleFoam" | tail -4
