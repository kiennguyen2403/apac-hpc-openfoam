#!/bin/bash -l

#PBS -N openfoam_clang_build
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=01:00:00
#PBS -l jobfs=40GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -j oe

set -o pipefail

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
OPENFOAM_DIR="$OPENFOAM_ROOT/OpenFOAM-v2512"
LOG_DIR="$OPENFOAM_ROOT/logs/build/v2512-clang-openmpi"
BUILD_JOBS="${PBS_NCPUS:-32}"
EXPECTED_COMMIT=87ed40d256d22ea38fcc648dfc82a22162427b18

mkdir -p "$LOG_DIR"
exec > >(tee "$LOG_DIR/build-${PBS_JOBID:-manual}.log") 2>&1

[[ $(git -C "$OPENFOAM_DIR" rev-parse HEAD) == "$EXPECTED_COMMIT" ]] || {
    echo "ERROR: OpenFOAM-v2512 source is not pinned at the expected commit" >&2
    echo "expected=$EXPECTED_COMMIT actual=$(git -C "$OPENFOAM_DIR" rev-parse HEAD)" >&2
    exit 1
}

module purge
module load clang/19.1.6
module load openmpi/5.0.8
module load metis/5.1.0
module load fftw3/3.3.10
module load boost/1.84.0

export CC=clang CXX=clang++ OMPI_CC=clang OMPI_CXX=clang++
export WM_COMPILER_TYPE=system WM_COMPILER=Clang WM_CC=clang WM_CXX=clang++
export WM_MPLIB=SYSTEMOPENMPI
# wmake consumes FOAM_EXTRA_*FLAGS, not the exported WM_CXXFLAGS.
# Match the AVX-512 instruction-set target of ICX's -xCORE-AVX512.
export FOAM_EXTRA_CFLAGS="-O3 -march=skylake-avx512"
export FOAM_EXTRA_CXXFLAGS="$FOAM_EXTRA_CFLAGS"

# Command-line settings are evaluated after etc/prefs.sh and therefore
# override this workspace's default ICX/Intel MPI preferences.
source "$OPENFOAM_DIR/etc/bashrc" \
    WM_COMPILER_TYPE=system WM_COMPILER=Clang WM_CC=clang WM_CXX=clang++ \
    WM_MPLIB=SYSTEMOPENMPI WM_COMPILE_OPTION=Opt
set -eE
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512

export CC=clang CXX=clang++ OMPI_CC=clang OMPI_CXX=clang++
export WM_MPLIB=SYSTEMOPENMPI FOAM_MPI=sys-openmpi MPI_ARCH_PATH="$OPENMPI_ROOT"
export FOAM_EXTRA_CFLAGS="-O3 -march=skylake-avx512"
export FOAM_EXTRA_CXXFLAGS="$FOAM_EXTRA_CFLAGS"
export WM_NCOMPPROCS="$BUILD_JOBS"
# Object/dependency files stay on the build node, leaving scratch inodes free.
export FOAM_BUILDROOT="${PBS_JOBFS:?Submit through PBS with jobfs}/clang-objects"
mkdir "$FOAM_BUILDROOT"
export FOAM_MODULE_PREFIX=none

[[ "$WM_OPTIONS" == linux64ClangDPInt32Opt ]] || {
    echo "ERROR: unexpected WM_OPTIONS: $WM_OPTIONS" >&2
    exit 1
}

export FFTW_ARCH_PATH="${FFTW3_ROOT:-/apps/fftw3/3.3.10}"
export FFTW_INC_DIR="$FFTW_ARCH_PATH/include"
export FFTW_LIB_DIR="$FFTW_ARCH_PATH/lib/GNU"
export LIBRARY_PATH="$FFTW_LIB_DIR:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$FFTW_LIB_DIR:${LD_LIBRARY_PATH:-}"

{
    echo "source_commit=$(git -C "$OPENFOAM_DIR" rev-parse HEAD)"
    echo "source_tag=OpenFOAM-v2512"
    echo "wm_project_version=$WM_PROJECT_VERSION"
    echo "wm_options=$WM_OPTIONS"
    echo "wm_mplib=$WM_MPLIB"
    echo "foam_mpi=$FOAM_MPI"
    echo "compiler=$($CXX --version | head -1)"
    echo "mpi=$(mpirun --version | head -1)"
    echo "build_jobs=$BUILD_JOBS"
    echo "cxxflags=$(wmake -show-cxxflags)"
    echo "foam_buildroot=$FOAM_BUILDROOT"
    echo "foam_extra_cxxflags=$FOAM_EXTRA_CXXFLAGS"
    echo "timestamp_utc=$(date -u +%FT%TZ)"
} | tee "$LOG_DIR/build-info-${PBS_JOBID:-manual}.txt"
module list > "$LOG_DIR/modules-${PBS_JOBID:-manual}.txt" 2>&1

cd "$WM_PROJECT_DIR"
# Dependencies are supplied by read-only site modules; do not rebuild ThirdParty.
wmake/src/Allmake
src/Allwmake -s -j "$BUILD_JOBS" -q
applications/Allwmake -s -j "$BUILD_JOBS" -q

for executable in simpleFoam decomposePar renumberMesh potentialFoam applyBoundaryLayer blockMesh foamDictionary; do
    executable_path=$(command -v "$executable")
    [[ $executable_path == "$OPENFOAM_DIR"/* ]]
    "$executable" -help >/dev/null 2>&1
    echo "$executable=$executable_path"
done

ldd "$(command -v simpleFoam)" > "$LOG_DIR/simpleFoam-ldd-${PBS_JOBID:-manual}.txt"
! grep -q 'not found' "$LOG_DIR/simpleFoam-ldd-${PBS_JOBID:-manual}.txt"
echo "BUILD_SUCCESS"
