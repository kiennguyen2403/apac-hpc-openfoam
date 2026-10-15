#!/bin/bash -l

#PBS -N openfoam_build
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=01:00:00
#PBS -l wd
#PBS -o /scratch/ya41/kn3896/openfoam/logs/build/openfoam_build.out
#PBS -e /scratch/ya41/kn3896/openfoam/logs/build/openfoam_build.err

set -o pipefail

# ============================================================
# Paths
# ============================================================

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
OPENFOAM_DIR="$OPENFOAM_ROOT/OpenFOAM-v2512"
THIRDPARTY_DIR="$OPENFOAM_ROOT/ThirdParty-v2512"
BUILD_LOG_DIR="$OPENFOAM_ROOT/build"
BUILD_JOBS="${PBS_NCPUS:-32}"

mkdir -p "$BUILD_LOG_DIR"
LOG_FILE="$BUILD_LOG_DIR/openfoam_build_${PBS_JOBID:-manual}.out"
exec > >(tee "$LOG_FILE") 2>&1

# Show the exact command/line if the script fails.
trap 'rc=$?; echo "ERROR: line $LINENO: $BASH_COMMAND (exit $rc)" >&2' ERR

# ============================================================
# Modules
# ============================================================

module purge
module load intel-compiler-llvm/2025.3.2
module load intel-mpi/2021.17.2
module load metis/5.1.0
module load intel-mkl/2025.3.1
module load boost/1.84.0

# ============================================================
# Intel compiler / MPI
# ============================================================

export CC=icx
export CXX=icpx

export WM_COMPILER_TYPE=system
export WM_COMPILER=Icx
export WM_CC=icx
export WM_CXX=icpx
export WM_MPLIB=INTELMPI

export I_MPI_CC=icx
export I_MPI_CXX=icpx

export FOAM_MPI=intelmpi
export MPI_ARCH_PATH="$I_MPI_ROOT"

export WM_CFLAGS="-O2 -fp-model=fast -m64 -fPIC"
export WM_CXXFLAGS="$WM_CFLAGS"

export CFLAGS="-O2 -fPIC"
export CXXFLAGS="$CFLAGS"

export WM_NCOMPPROCS="$BUILD_JOBS"

# ============================================================
# Load OpenFOAM
# ============================================================

echo "Loading OpenFOAM..."
source "$OPENFOAM_DIR/etc/bashrc" || true
set -eE
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR
export WM_THIRD_PARTY_DIR="$THIRDPARTY_DIR"

if [[ -z "${WM_PROJECT_DIR:-}" ]]; then
    echo "ERROR: OpenFOAM environment was not initialized"
    exit 1
fi

# OpenFOAM bashrc may overwrite compiler/MPI settings, so set them again.
export CC=icx
export CXX=icpx

export WM_COMPILER_TYPE=system
export WM_COMPILER=Icx
export WM_CC=icx
export WM_CXX=icpx
export WM_MPLIB=INTELMPI

export I_MPI_CC=icx
export I_MPI_CXX=icpx

export FOAM_MPI=intelmpi
export MPI_ARCH_PATH="$I_MPI_ROOT"

export WM_CFLAGS="-O2 -fp-model=fast -m64 -fPIC"
export WM_CXXFLAGS="$WM_CFLAGS"

export CFLAGS="-O2 -fPIC"
export CXXFLAGS="$CFLAGS"

export WM_NCOMPPROCS="$BUILD_JOBS"

# ============================================================
# PETSc + Intel MKL configuration
# ============================================================

PETSC_INSTALL_DIR=""
PETSC_LIB_DIR=""

for cand in \
    "$THIRDPARTY_DIR/platforms/linux64IcxDPInt32/petsc-3.21.2" \
    "$THIRDPARTY_DIR/platforms/linux64IcxDPInt32/petsc" \
    "$THIRDPARTY_DIR/petsc-3.21.2" \
    "$THIRDPARTY_DIR/petsc"
 do
    if [[ -f "$cand/lib/libpetsc.so" ]]; then
        PETSC_INSTALL_DIR="$cand"
        PETSC_LIB_DIR="$cand/lib"
        break
    fi
done

if [[ -z "$PETSC_INSTALL_DIR" || -z "$PETSC_LIB_DIR" ]]; then
    echo "ERROR: PETSc library not found in any expected location."
    echo "Expected candidates:"
    for cand in \
        "$THIRDPARTY_DIR/platforms/linux64IcxDPInt32/petsc-3.21.2" \
        "$THIRDPARTY_DIR/platforms/linux64IcxDPInt32/petsc" \
        "$THIRDPARTY_DIR/petsc-3.21.2" \
        "$THIRDPARTY_DIR/petsc"
    do
        echo "  $cand/lib/libpetsc.so"
    done
    echo "       Build the PETSc+MKL stack first using scripts/build/math_libraries/build_petsc_mkl.sh"
    exit 1
fi

export PETSC_DIR="$PETSC_INSTALL_DIR"
export PETSC_ARCH=""

export LD_LIBRARY_PATH="$PETSC_LIB_DIR:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="$PETSC_LIB_DIR:${LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$PETSC_LIB_DIR/pkgconfig:${PKG_CONFIG_PATH:-}"

echo
echo "============================================================"
echo "Using PETSc + Intel MKL"
echo "============================================================"
echo "PETSC_DIR = $PETSC_DIR"
echo "PETSC_ARCH = $PETSC_ARCH"
echo "PETSc lib = $PETSC_LIB_DIR/libpetsc.so"

echo
echo "WM_PROJECT_DIR     = $WM_PROJECT_DIR"
echo "WM_THIRD_PARTY_DIR = $WM_THIRD_PARTY_DIR"
echo "WM_OPTIONS        = $WM_OPTIONS"

# ============================================================
# FFTW3 interface backed by Intel MKL
# ============================================================

export FFTW_ARCH_PATH="$THIRDPARTY_DIR/platforms/linux64IcxDPInt32/fftw-mkl-2025.3.1"
export FFTW_INC_DIR="$FFTW_ARCH_PATH/include"
export FFTW_LIB_DIR="$FFTW_ARCH_PATH/lib"

if [[ ! -f "$FFTW_INC_DIR/fftw3.h" || ! -f "$FFTW_LIB_DIR/libfftw3.so" ]]; then
    echo "ERROR: Intel MKL FFTW3 interface was not found under:"
    echo "       $FFTW_ARCH_PATH"
    echo "       Build it first using build_fftw_mkl.sh"
    exit 1
fi

export CPATH="$FFTW_INC_DIR:${CPATH:-}"
export LIBRARY_PATH="$FFTW_LIB_DIR:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$FFTW_LIB_DIR:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$FFTW_LIB_DIR/pkgconfig:${PKG_CONFIG_PATH:-}"

echo
echo "============================================================"
echo "Using FFTW3 interface backed by Intel MKL"
echo "============================================================"
echo "FFTW_ARCH_PATH = $FFTW_ARCH_PATH"
echo "FFTW include   = $FFTW_INC_DIR/fftw3.h"
echo "FFTW library   = $FFTW_LIB_DIR/libfftw3.so"

# ============================================================
# Build OpenFOAM
# ============================================================

cd "$WM_PROJECT_DIR"

echo
echo "=========================================="
echo "Building OpenFOAM with $BUILD_JOBS cores"
echo "=========================================="

./Allwmake \
    -prefix=openfoam \
    -s \
    -j "$BUILD_JOBS" \
    -q

# ============================================================
# Basic verification
# ============================================================

echo
echo "=========================================="
echo "Build finished"
echo "=========================================="

echo "WM_OPTIONS = $WM_OPTIONS"

echo
echo "OpenFOAM executables:"
command -v blockMesh || true
command -v decomposePar || true
command -v foamRun || true

PETSC_FOAM_LIB="$FOAM_LIBBIN/libpetscFoam.so"

if [[ ! -f "$PETSC_FOAM_LIB" ]]; then
    echo "ERROR: PETSc/OpenFOAM solver module was not built: $PETSC_FOAM_LIB"
    echo "       Ensure modules/external-solver is initialized"
    exit 1
fi

PETSC_FOAM_LDD="$(ldd "$PETSC_FOAM_LIB")"
! grep -q 'not found' <<< "$PETSC_FOAM_LDD" || {
    echo "$PETSC_FOAM_LDD"
    exit 1
}
grep -q 'libpetsc' <<< "$PETSC_FOAM_LDD" || {
    echo "ERROR: libpetscFoam is not linked to PETSc"
    exit 1
}
grep -q 'libmkl_core' <<< "$PETSC_FOAM_LDD" || {
    echo "ERROR: libpetscFoam/PETSc is not linked to Intel MKL"
    exit 1
}

echo
echo "PETSc/OpenFOAM module: $PETSC_FOAM_LIB"

echo
echo "=========================================="
echo "OpenFOAM build completed successfully"
echo "=========================================="
