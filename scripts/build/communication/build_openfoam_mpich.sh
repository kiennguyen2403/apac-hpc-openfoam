#!/bin/bash -l
#PBS -N of2512_icx_mpich
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/ya41
#PBS -l jobfs=20GB
#PBS -l wd
#PBS -j oe

set -o pipefail

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
OPENFOAM_DIR=$OPENFOAM_ROOT/OpenFOAM-v2512
THIRDPARTY_DIR=$OPENFOAM_ROOT/ThirdParty-v2512
MPICH_VERSION=mpich-3.3
MPICH_ARCHIVE=$THIRDPARTY_DIR/downloads/$MPICH_VERSION.tar.gz
BUILD_JOBS=${PBS_NCPUS:-32}

module purge
module load intel-compiler-llvm/2025.3.2
module load metis/5.1.0
module load fftw3/3.3.10
module load boost/1.84.0

export CC=icx CXX=icpx
source "$OPENFOAM_DIR/etc/bashrc" \
    WM_COMPILER_TYPE=system \
    WM_COMPILER=Icx \
    WM_CC=icx \
    WM_CXX=icpx \
    WM_MPLIB=MPICH \
    WM_COMPILE_OPTION=OptAVX512

set -eE
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR

[[ $WM_PROJECT_VERSION == v2512 ]]
[[ $WM_OPTIONS == linux64IcxDPInt32OptAVX512 ]]
[[ -s $MPICH_ARCHIVE ]]
: "${PBS_JOBFS:?This build must run inside a PBS job}"

# Source trees are unpacked into node-local jobfs to avoid scratch inode quota.
# Only installed MPICH files are retained in the writable v2512 platforms tree.
RUNTIME_TP=$PBS_JOBFS/ThirdParty-v2512
mkdir -p "$RUNTIME_TP/sources/mpich/$MPICH_VERSION"
cp "$THIRDPARTY_DIR/makeMPICH" "$RUNTIME_TP/"
cp -a "$THIRDPARTY_DIR/etc" "$RUNTIME_TP/"
ln -s "$THIRDPARTY_DIR/platforms" "$RUNTIME_TP/platforms"
tar -xzf "$MPICH_ARCHIVE" --strip-components=1 \
    -C "$RUNTIME_TP/sources/mpich/$MPICH_VERSION"

export WM_THIRD_PARTY_DIR=$RUNTIME_TP
export WM_NCOMPPROCS=$BUILD_JOBS
export CC=icx CXX=icpx
export WM_MPLIB=MPICH
export FOAM_MPI=$MPICH_VERSION
unset F77 FFLAGS F90 F90FLAGS FC FCFLAGS

cd "$WM_THIRD_PARTY_DIR"
./makeMPICH -force "$MPICH_VERSION" -- --disable-opencl --with-device=ch3:nemesis

export WM_THIRD_PARTY_DIR=$THIRDPARTY_DIR
export MPI_ARCH_PATH=$THIRDPARTY_DIR/platforms/$WM_ARCH$WM_COMPILER/$MPICH_VERSION
export MPI_HOME=$MPI_ARCH_PATH
export PATH=$MPI_ARCH_PATH/bin:$PATH
export LD_LIBRARY_PATH=$MPI_ARCH_PATH/lib:${LD_LIBRARY_PATH:-}
export LIBRARY_PATH=$MPI_ARCH_PATH/lib:${LIBRARY_PATH:-}
export PTSCOTCH_ARCH_PATH=none

[[ -x $MPI_ARCH_PATH/bin/mpiexec ]]
"$MPI_ARCH_PATH/bin/mpiexec" --version

cd "$OPENFOAM_DIR"
# Only Pstream is MPI-specific for the existing v2512 application build.
# Avoid top-level Allwmake: it also rebuilds optional ThirdParty packages.
(cd src/Pstream && ./Allwmake -j "$BUILD_JOBS")
MPI_PSTREAM_LIB=$FOAM_LIBBIN/$FOAM_MPI/libPstream.so
[[ -s $MPI_PSTREAM_LIB ]]

for executable in simpleFoam decomposePar; do
    executable_path=$(command -v "$executable")
    [[ $executable_path == "$OPENFOAM_DIR"/* ]]
    "$executable" -help >/dev/null 2>&1
done
echo BUILD_SUCCESS
