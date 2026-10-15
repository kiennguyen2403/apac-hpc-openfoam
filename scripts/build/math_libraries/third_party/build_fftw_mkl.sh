#!/bin/bash -l

#PBS -N fftw_mkl_build
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=00:30:00
#PBS -l wd
#PBS -j oe
#PBS -o /scratch/ya41/kn3896/openfoam/logs/fftw_mkl/mkl_build.out
#PBS -e /scratch/ya41/kn3896/openfoam/logs/fftw_mkl/mkl_build.err

set -eEo pipefail

# Build Intel oneMKL's FFTW3 C compatibility interface.  This is not a build
# of the upstream FFTW implementation: applications call the FFTW3 API and
# the wrapper dispatches the transforms to oneMKL.

# ============================================================
# Paths
# ============================================================

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
INSTALL_DIR="$OPENFOAM_ROOT/ThirdParty-v2512/platforms/linux64IcxDPInt32/fftw-mkl-2025.3.1"
BUILD_ROOT="$OPENFOAM_ROOT/build/fftw-mkl"
BUILD_DIR="$BUILD_ROOT/${PBS_JOBID:-manual}"
BUILD_LOG_DIR="$OPENFOAM_ROOT/build"
BUILD_JOBS="${PBS_NCPUS:-32}"

mkdir -p "$INSTALL_DIR/include" "$INSTALL_DIR/lib/pkgconfig" "$BUILD_DIR" "$BUILD_LOG_DIR"
LOG_FILE="$BUILD_LOG_DIR/fftw_mkl_build_${PBS_JOBID:-manual}.out"
exec > >(tee "$LOG_FILE") 2>&1

trap 'rc=$?; echo "ERROR: line $LINENO: $BASH_COMMAND (exit $rc)" >&2' ERR

# ============================================================
# Modules
# ============================================================

module purge
module load intel-compiler-llvm/2025.3.2
module load intel-mkl/2025.3.1

export CC=icx

if [[ -z "${MKLROOT:-}" || ! -d "$MKLROOT" ]]; then
    echo "ERROR: MKLROOT is not set to an existing directory"
    exit 1
fi

WRAPPER_DIR="$MKLROOT/share/mkl/interfaces/fftw3xc"
WRAPPER_MAKEFILE="$WRAPPER_DIR/makefile"

if [[ ! -f "$WRAPPER_MAKEFILE" ]]; then
    echo "ERROR: oneMKL FFTW3 C wrapper makefile was not found:"
    echo "       $WRAPPER_MAKEFILE"
    exit 1
fi

echo "============================================================"
echo "Building the FFTW3 interface to Intel oneMKL"
echo "============================================================"
echo "MKLROOT     = $MKLROOT"
echo "INSTALL_DIR = $INSTALL_DIR"
echo "BUILD_DIR   = $BUILD_DIR"
echo "BUILD_JOBS  = $BUILD_JOBS"

# ============================================================
# Build the Intel FFTW3 wrapper archive
# ============================================================

# Intel's makefile removes its temporary object directory after archiving.
# INSTALL_DIR is therefore set to the per-job build directory first.  -fPIC
# allows the wrapper archive to be incorporated into libfftw3.so below.
make -C "$WRAPPER_DIR" \
    -f "$WRAPPER_MAKEFILE" \
    -j "$BUILD_JOBS" \
    libintel64 \
    compiler=intel \
    MKLROOT="$MKLROOT" \
    INSTALL_DIR="$BUILD_DIR" \
    INSTALL_LIBNAME=libfftw3xc_intel.a \
    CFLAGS="-O2 -fPIC"

test -s "$BUILD_DIR/libfftw3xc_intel.a"

# ============================================================
# Install FFTW-compatible headers and libraries
# ============================================================

install -m 0644 "$MKLROOT/include/fftw/fftw3.h" "$INSTALL_DIR/include/fftw3.h"
install -m 0644 "$MKLROOT/include/fftw/fftw3_mkl.h" "$INSTALL_DIR/include/fftw3_mkl.h"
install -m 0644 "$BUILD_DIR/libfftw3xc_intel.a" "$INSTALL_DIR/lib/libfftw3xc_intel.a"

# OpenFOAM normally links FFTW as -lfftw3.  Give the MKL wrapper that standard
# library name and retain all wrapper entry points when creating the DSO.
icx -shared \
    -Wl,-soname,libfftw3.so \
    -Wl,--whole-archive "$BUILD_DIR/libfftw3xc_intel.a" -Wl,--no-whole-archive \
    -L"$MKLROOT/lib" -Wl,-rpath,"$MKLROOT/lib" \
    -Wl,--no-as-needed -lmkl_rt \
    -lpthread -lm -ldl \
    -o "$INSTALL_DIR/lib/libfftw3.so"

# The wrapper archive contains both fftw_* (double) and fftwf_* (float)
# symbols.  These links support software that requests the conventional names.
ln -sfn libfftw3xc_intel.a "$INSTALL_DIR/lib/libfftw3.a"
ln -sfn libfftw3.so "$INSTALL_DIR/lib/libfftw3f.so"
ln -sfn libfftw3xc_intel.a "$INSTALL_DIR/lib/libfftw3f.a"

cat > "$BUILD_DIR/fftw3.pc" <<EOF
prefix=$INSTALL_DIR
exec_prefix=\${prefix}
libdir=\${prefix}/lib
includedir=\${prefix}/include

Name: fftw3-mkl
Description: FFTW3 C compatibility interface backed by Intel oneMKL
Version: 3.x-mkl-2025.3.1
Libs: -L\${libdir} -Wl,-rpath,\${libdir} -lfftw3
Libs.private: -L$MKLROOT/lib -lmkl_rt -lpthread -lm -ldl
Cflags: -I\${includedir}
EOF
install -m 0644 "$BUILD_DIR/fftw3.pc" "$INSTALL_DIR/lib/pkgconfig/fftw3.pc"

# ============================================================
# Numerical smoke test
# ============================================================

cat > "$BUILD_DIR/test_fftw_mkl.c" <<'EOF'
#include <math.h>
#include <stdio.h>
#include <fftw3.h>

int main(void)
{
    double in[4] = {1.0, 2.0, 3.0, 4.0};
    fftw_complex out[3];
    fftw_plan plan = fftw_plan_dft_r2c_1d(4, in, out, FFTW_ESTIMATE);

    if (!plan) return 2;
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    if (fabs(out[0][0] - 10.0) > 1e-12 || fabs(out[0][1]) > 1e-12) {
        fprintf(stderr, "unexpected DC result: %.17g %+.17gi\n", out[0][0], out[0][1]);
        return 3;
    }

    printf("FFTW3/MKL smoke test passed: DC = %.1f %+.1fi\n", out[0][0], out[0][1]);
    return 0;
}
EOF

icx -O2 \
    -I"$INSTALL_DIR/include" \
    "$BUILD_DIR/test_fftw_mkl.c" \
    -L"$INSTALL_DIR/lib" -Wl,-rpath,"$INSTALL_DIR/lib" \
    -lfftw3 -o "$BUILD_DIR/test_fftw_mkl"

"$BUILD_DIR/test_fftw_mkl"

echo
echo "MKL dependency:"
ldd "$INSTALL_DIR/lib/libfftw3.so" | grep -E 'mkl|not found' || true

echo
echo "============================================================"
echo "FFTW3 interface to Intel oneMKL built successfully"
echo "============================================================"
echo "FFTW_ARCH_PATH=$INSTALL_DIR"
echo "FFTW_INC_DIR=$INSTALL_DIR/include"
echo "FFTW_LIB_DIR=$INSTALL_DIR/lib"
echo "PKG_CONFIG_PATH=$INSTALL_DIR/lib/pkgconfig:\${PKG_CONFIG_PATH:-}"
echo "LD_LIBRARY_PATH=$INSTALL_DIR/lib:\${LD_LIBRARY_PATH:-}"

