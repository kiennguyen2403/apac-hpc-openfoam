#!/bin/bash -l

#PBS -N openfoam_openmpi_build
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -o /scratch/ya41/kn3896/openfoam/logs/build/openfoam_openmpi_build.out
#PBS -e /scratch/ya41/kn3896/openfoam/logs/build/openfoam_openmpi_build.err

set -o pipefail

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
OPENFOAM_DIR="$OPENFOAM_ROOT/OpenFOAM-v2512"
BUILD_JOBS="${PBS_NCPUS:-32}"

module purge
module load intel-compiler-llvm/2025.3.2
module load openmpi/5.0.8
module load metis/5.1.0
module load fftw3/3.3.10
module load intel-mkl/2025.3.1
module load boost/1.84.0

# Keep the compiler identical to the Intel MPI build so only MPI changes.
export CC=icx CXX=icpx
export OMPI_CC=icx OMPI_CXX=icpx
export WM_COMPILER_TYPE=system WM_COMPILER=Icx
export WM_CC=icx WM_CXX=icpx
export WM_MPLIB=SYSTEMOPENMPI WM_COMPILE_OPTION=OptAVX512
export WM_CFLAGS="-O2 -fp-model=fast -m64 -fPIC"
export WM_CXXFLAGS="$WM_CFLAGS"
export CFLAGS="-O2 -fPIC" CXXFLAGS="-O2 -fPIC"

source "$OPENFOAM_DIR/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx WM_MPLIB=SYSTEMOPENMPI WM_COMPILE_OPTION=OptAVX512
set -eE
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512
[[ $WM_OPTIONS == linux64IcxDPInt32OptAVX512 ]] || { echo "ERROR: unexpected WM_OPTIONS=$WM_OPTIONS" >&2; exit 2; }
[[ -n "${WM_PROJECT_DIR:-}" ]] || {
    echo "ERROR: OpenFOAM environment was not initialized" >&2
    exit 1
}

# bashrc contains defaults, so enforce this build's settings afterwards.
export CC=icx CXX=icpx
export OMPI_CC=icx OMPI_CXX=icpx
export WM_COMPILER_TYPE=system WM_COMPILER=Icx
export WM_CC=icx WM_CXX=icpx
export WM_MPLIB=SYSTEMOPENMPI WM_COMPILE_OPTION=OptAVX512
export FOAM_MPI=sys-openmpi
export MPI_ARCH_PATH="$OPENMPI_ROOT"
export WM_CFLAGS="-O2 -fp-model=fast -m64 -fPIC"
export WM_CXXFLAGS="$WM_CFLAGS"
export CFLAGS="-O2 -fPIC" CXXFLAGS="-O2 -fPIC"
export WM_NCOMPPROCS="$BUILD_JOBS"

export FFTW_ARCH_PATH="${FFTW3_ROOT:-/apps/fftw3/3.3.10}"
export FFTW_INC_DIR="$FFTW_ARCH_PATH/include"
export FFTW_LIB_DIR="$FFTW_ARCH_PATH/lib/Intel"
export LIBRARY_PATH="$FFTW_LIB_DIR:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$FFTW_LIB_DIR:${LD_LIBRARY_PATH:-}"

echo "Compiler: $($CXX --version | head -1)"
echo "MPI:      $(mpirun --version | head -1)"
echo "WM_OPTIONS=$WM_OPTIONS, WM_MPLIB=$WM_MPLIB, FOAM_MPI=$FOAM_MPI"


cd "$WM_PROJECT_DIR"
./Allwmake -s -j "$BUILD_JOBS" -q

for command_name in blockMesh decomposePar foamRun mpirun; do
    command -v "$command_name" || true
done

echo "OpenFOAM OpenMPI build completed successfully"

