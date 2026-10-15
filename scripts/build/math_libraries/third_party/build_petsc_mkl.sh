#!/bin/bash -l
#PBS -N of2512_petsc_mkl
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/ya41
#PBS -l jobfs=30GB
#PBS -l wd
#PBS -j oe

set -o pipefail

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
OPENFOAM_DIR=$OPENFOAM_ROOT/OpenFOAM-v2512
THIRDPARTY_DIR=$OPENFOAM_ROOT/ThirdParty-v2512
PETSC_VERSION=petsc-3.21.2
HYPRE_VERSION=hypre-2.31.0
PETSC_ARCHIVE=$THIRDPARTY_DIR/downloads/$PETSC_VERSION.tar.gz
HYPRE_ARCHIVE=$THIRDPARTY_DIR/downloads/$HYPRE_VERSION.tar.gz
BUILD_JOBS=${PBS_NCPUS:-32}

module purge
module load intel-compiler-llvm/2025.3.2
module load intel-mpi/2021.17.2
module load metis/5.1.0
module load intel-mkl/2025.3.1
module load boost/1.84.0

export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
source "$OPENFOAM_DIR/etc/bashrc" WM_COMPILER=Icx WM_MPLIB=INTELMPI

set -eE
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR

[[ $WM_PROJECT_VERSION == v2512 ]]
for archive in "$PETSC_ARCHIVE" "$HYPRE_ARCHIVE"; do
    [[ -s $archive ]]
done
: "${PBS_JOBFS:?This build must run inside a PBS job}"

# Unpack large sources into node-local jobfs. Build/install products pass
# through the platforms symlink and remain in the writable scratch tree.
RUNTIME_TP=$PBS_JOBFS/ThirdParty-v2512
mkdir -p \
    "$RUNTIME_TP/sources/petsc/$PETSC_VERSION" \
    "$RUNTIME_TP/sources/hypre/$HYPRE_VERSION"
cp "$THIRDPARTY_DIR/makePETSC" "$THIRDPARTY_DIR/makeHYPRE" "$RUNTIME_TP/"
cp -a "$THIRDPARTY_DIR/etc" "$RUNTIME_TP/"
ln -s "$THIRDPARTY_DIR/platforms" "$RUNTIME_TP/platforms"
tar -xzf "$PETSC_ARCHIVE" --strip-components=1 \
    -C "$RUNTIME_TP/sources/petsc/$PETSC_VERSION"
tar -xzf "$HYPRE_ARCHIVE" --strip-components=1 \
    -C "$RUNTIME_TP/sources/hypre/$HYPRE_VERSION"

export WM_THIRD_PARTY_DIR=$RUNTIME_TP
export WM_NCOMPPROCS=$BUILD_JOBS
export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_THREADING_LAYER=SEQUENTIAL
export LIBRARY_PATH=$MKLROOT/lib:${LIBRARY_PATH:-}
export LD_LIBRARY_PATH=$MKLROOT/lib:${LD_LIBRARY_PATH:-}

cd "$WM_THIRD_PARTY_DIR"
./makeHYPRE -force "$HYPRE_VERSION"
./makePETSC -force "$HYPRE_VERSION" "$PETSC_VERSION" -- \
    --with-blaslapack-dir="$MKLROOT"

export WM_THIRD_PARTY_DIR=$THIRDPARTY_DIR
PETSC_LIB=$(find "$THIRDPARTY_DIR/platforms" \
    \( -type f -o -type l \) -name libpetsc.so -print -quit)
[[ -n $PETSC_LIB ]]
if ldd "$PETSC_LIB" | grep -q 'not found'; then
    ldd "$PETSC_LIB" | grep 'not found' >&2
    exit 1
fi
echo "PETSc library: $PETSC_LIB"
echo BUILD_SUCCESS
