#!/bin/bash -l

#PBS -N openfoam_gcc_build
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=01:00:00
#PBS -l wd
#PBS -o /scratch/ya41/kn3896/openfoam/logs/build/openfoam_gcc_build.out
#PBS -e /scratch/ya41/kn3896/openfoam/logs/build/openfoam_gcc_build.err

set -o pipefail

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
OPENFOAM_DIR="$OPENFOAM_ROOT/OpenFOAM-v2512"
BUILD_JOBS="${PBS_NCPUS:-32}"

module purge
module load gcc/14.2.0
module load openmpi/5.0.8
module load metis/5.1.0
module load fftw3/3.3.10
module load boost/1.84.0

export CC=gcc CXX=g++ OMPI_CC=gcc OMPI_CXX=g++
export WM_COMPILER_TYPE=system WM_COMPILER=Gcc WM_CC=gcc WM_CXX=g++
export WM_MPLIB=SYSTEMOPENMPI
export WM_CFLAGS="-O2 -m64 -fPIC"
export WM_CXXFLAGS="$WM_CFLAGS"
export CFLAGS="-O2 -m64 -fPIC" CXXFLAGS="-O2 -m64 -fPIC"

# Command-line settings are evaluated after etc/prefs.sh and therefore
# override this workspace's default ICX/Intel MPI preferences.
source "$OPENFOAM_DIR/etc/bashrc" \
    WM_COMPILER_TYPE=system WM_COMPILER=Gcc WM_CC=gcc WM_CXX=g++ \
    WM_MPLIB=SYSTEMOPENMPI
set -eE
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512

export CC=gcc CXX=g++ OMPI_CC=gcc OMPI_CXX=g++
export WM_MPLIB=SYSTEMOPENMPI FOAM_MPI=sys-openmpi MPI_ARCH_PATH="$OPENMPI_ROOT"
export WM_CFLAGS="-O2 -m64 -fPIC"
export WM_CXXFLAGS="$WM_CFLAGS"
export CFLAGS="-O2 -m64 -fPIC" CXXFLAGS="-O2 -m64 -fPIC"
export WM_NCOMPPROCS="$BUILD_JOBS"

[[ "$WM_OPTIONS" == linux64GccDPInt32Opt ]] || {
    echo "ERROR: unexpected WM_OPTIONS: $WM_OPTIONS" >&2
    exit 1
}

export FFTW_ARCH_PATH="${FFTW3_ROOT:-/apps/fftw3/3.3.10}"
export FFTW_INC_DIR="$FFTW_ARCH_PATH/include"
export FFTW_LIB_DIR="$FFTW_ARCH_PATH/lib/GNU"
export LIBRARY_PATH="$FFTW_LIB_DIR:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$FFTW_LIB_DIR:${LD_LIBRARY_PATH:-}"

echo "Compiler: $($CXX --version | head -1)"
echo "MPI:      $(mpirun --version | head -1)"
echo "WM_OPTIONS=$WM_OPTIONS, WM_MPLIB=$WM_MPLIB, FOAM_MPI=$FOAM_MPI"


cd "$WM_PROJECT_DIR"
./Allwmake -s -j "$BUILD_JOBS" -q

for command_name in blockMesh decomposePar foamRun simpleFoam mpirun; do
    command -v "$command_name" || true
done

echo "OpenFOAM GCC/OpenMPI build completed successfully"

