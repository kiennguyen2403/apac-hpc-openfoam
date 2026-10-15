#!/bin/bash -l
#PBS -N of2512_icx_impi
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=96GB
#PBS -l walltime=02:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -j oe

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
OPENFOAM_DIR=$OPENFOAM_ROOT/OpenFOAM-v2512
LOG_DIR=$OPENFOAM_ROOT/logs/build/v2512-icx-intelmpi
BUILD_JOBS=${PBS_NCPUS:-48}
EXPECTED_COMMIT=87ed40d256d22ea38fcc648dfc82a22162427b18

mkdir -p "$LOG_DIR"
exec > >(tee "$LOG_DIR/build-${PBS_JOBID:-manual}.log") 2>&1

module purge
module load intel-compiler-llvm/2025.3.2
module load intel-mpi/2021.17.2
module load metis/5.1.0
module load fftw3/3.3.10
module load boost/1.84.0

export CC=icx
export CXX=icpx
export I_MPI_CC=icx
export I_MPI_CXX=icpx

source "$OPENFOAM_DIR/etc/bashrc" \
    WM_COMPILER=Icx \
    WM_MPLIB=INTELMPI \
    WM_COMPILE_OPTION=OptAVX512

# The OpenFOAM environment loader probes optional files with nonzero return
# codes. Enable strict error handling only after it has finished, and never
# enable nounset (-u) for OpenFOAM shell environments.
set -eEo pipefail
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR

# Use the site's read-only v2512 dependency tree. Serial Scotch has a C ABI
# and can be used by ICX. Its PT-Scotch is OpenMPI-linked, so disable it for
# this Intel MPI build.
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512
export PTSCOTCH_ARCH_PATH=none
export WM_NCOMPPROCS=$BUILD_JOBS
export I_MPI_LIBRARY_KIND=release

[[ $WM_PROJECT_VERSION == v2512 ]]
[[ $WM_PRECISION_OPTION == DP ]]
[[ $WM_LABEL_SIZE == 32 ]]
[[ $WM_OPTIONS == linux64IcxDPInt32OptAVX512 ]]
[[ $WM_MPLIB == INTELMPI ]]
[[ $(git -C "$OPENFOAM_DIR" rev-parse HEAD) == "$EXPECTED_COMMIT" ]]

. "$WM_PROJECT_DIR/wmake/scripts/have_scotch"
have_scotch
if have_ptscotch; then
    echo "ERROR: OpenMPI PT-Scotch must not enter the Intel MPI build" >&2
    exit 2
fi

{
    echo "source_commit=$(git -C "$OPENFOAM_DIR" rev-parse HEAD)"
    echo "source_tag=OpenFOAM-v2512"
    echo "wm_project_version=$WM_PROJECT_VERSION"
    echo "wm_options=$WM_OPTIONS"
    echo "wm_mplib=$WM_MPLIB"
    echo "foam_mpi=$FOAM_MPI"
    echo "compiler=$(icpx --version | head -1)"
    echo "mpi=$(mpiexec --version 2>&1 | head -1)"
    echo "scotch_include=$SCOTCH_INC_DIR"
    echo "scotch_library=$SCOTCH_LIB_DIR"
    echo "build_jobs=$BUILD_JOBS"
    echo "timestamp_utc=$(date -u +%FT%TZ)"
} | tee "$LOG_DIR/build-info-${PBS_JOBID:-manual}.txt"
module list > "$LOG_DIR/modules-${PBS_JOBID:-manual}.txt" 2>&1

cd "$OPENFOAM_DIR"
./Allwmake -j "$BUILD_JOBS" -s -q

for executable in simpleFoam decomposePar renumberMesh potentialFoam applyBoundaryLayer; do
    executable_path=$(command -v "$executable")
    [[ $executable_path == "$OPENFOAM_DIR"/* ]]
    "$executable" -help >/dev/null 2>&1
    echo "$executable=$executable_path"
done

ldd "$(command -v simpleFoam)" > "$LOG_DIR/simpleFoam-ldd-${PBS_JOBID:-manual}.txt"
ldd "$(command -v decomposePar)" > "$LOG_DIR/decomposePar-ldd-${PBS_JOBID:-manual}.txt"
echo "BUILD_SUCCESS"
