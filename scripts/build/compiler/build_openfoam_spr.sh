#!/bin/bash -l
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=32,mem=64gb
#PBS -l walltime=01:40:00
#PBS -l jobfs=40GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -j oe

# Sapphire Rapids OpenFOAM-v2512 build for the one-node compiler comparison.
# Modelled on /scratch/ya41/apac-2026/openfoam/build/build.pbs:
# Intel MPI 2021.17.2, -O3 plus ISA flags via FOAM_EXTRA_*FLAGS, Allwmake -j 32.
# Submit once per compiler:
#   qsub -N of2512_spr_gcc   -v SPR_COMPILER=gcc   build_openfoam_spr.sh
#   qsub -N of2512_spr_clang -v SPR_COMPILER=clang build_openfoam_spr.sh
#   qsub -N of2512_spr_icx   -v SPR_COMPILER=icx   build_openfoam_spr.sh
# Each produces platforms/linux64<Compiler>DPInt32OptSPR (rules: wmake/rules/linux64*/c*OptSPR).
# Built on normalsr so the SPR executables can be smoke-tested on the build node.

set -o pipefail

OPENFOAM_ROOT=/scratch/ya41/kn3896/openfoam
OPENFOAM_DIR=$OPENFOAM_ROOT/OpenFOAM-v2512
BUILD_JOBS=${PBS_NCPUS:-32}
EXPECTED_COMMIT=87ed40d256d22ea38fcc648dfc82a22162427b18
: "${SPR_COMPILER:?set SPR_COMPILER=gcc|clang|icx}"
# Phase 1 flag tuning (icx only): SPR_TUNE=zmm|O2|ipo -> WM_COMPILE_OPTION=OptSPR<tune>
# (base -O level in wmake/rules/linux64Icx/c*OptSPR<tune>; extra flags below).
SPR_TUNE=${SPR_TUNE:-}
# MPI: SPR_MPI=intelmpi (default, Intel MPI 2021.17.2) or openmpi (Open MPI 4.1.7, the version used by
# module OpenFOAM/2512) -> WM_COMPILE_OPTION=OptSPR<tune>Ompi, so Intel MPI builds are not overwritten.
SPR_MPI=${SPR_MPI:-intelmpi}
case $SPR_MPI in
    intelmpi) MPI_MODULE=intel-mpi/2021.17.2; WM_MPLIB_SEL=INTELMPI;      OPT_SUFFIX= ;;
    openmpi)  MPI_MODULE=openmpi/4.1.7;       WM_MPLIB_SEL=SYSTEMOPENMPI; OPT_SUFFIX=Ompi ;;
    *) echo "ERROR: unknown SPR_MPI=$SPR_MPI" >&2; exit 2 ;;
esac
[[ -z $SPR_TUNE || $SPR_COMPILER == icx ]] || { echo "ERROR: SPR_TUNE is icx-only" >&2; exit 2; }
LOG_DIR=$OPENFOAM_ROOT/logs/build/v2512-$SPR_COMPILER-spr${SPR_TUNE:+-$SPR_TUNE}-$SPR_MPI

mkdir -p "$LOG_DIR"
exec > >(tee "$LOG_DIR/build-${PBS_JOBID:-manual}.log") 2>&1

module purge
case "$SPR_COMPILER" in
    gcc)
        # System as (binutils 2.30) rejects SPR's AVX512-FP16 instructions (e.g. vmovw).
        module load gcc/14.2.0 binutils/2.44
        WM_COMPILER=Gcc CC=gcc CXX=g++
        SPR_FLAGS="-march=sapphirerapids"
        ;;
    clang)
        module load clang/19.1.6
        WM_COMPILER=Clang CC=clang CXX=clang++
        SPR_FLAGS="-march=sapphirerapids -ffp-model=precise"
        ;;
    icx)
        module load intel-compiler-llvm/2025.3.2
        WM_COMPILER=Icx CC=icx CXX=icpx
        SPR_FLAGS="-xSAPPHIRERAPIDS -fp-model=precise"
        case "$SPR_TUNE" in
            '') ;;
            # Prefer 512-bit vectors. Masked zmm loops computed on inactive lanes and tripped
            # FOAM_SIGFPE in makeCellCentresAndVols (job 180741885); -fp-speculation=safe stops that.
            zmm) SPR_FLAGS+=" -qopt-zmm-usage=high -fp-speculation=safe" ;;
            # GCC-like multiply-add fusion across statements, other -fp-model=precise guarantees kept
            fma) SPR_FLAGS+=" -ffp-contract=fast" ;;
            O2)  ;;                                       # -O2 comes from the c*OptSPRO2 rule
            # Fat objects keep machine code so OpenFOAM's ld -r / ar steps still work;
            # wmake passes c++FLAGS to the link line, so -ipo also reaches every .so/exe link.
            # System ld 2.30's gold plugin cannot read ICX 2025 LTO bitcode ("Invalid encoding",
            # job 180741835): link with ICX's own lld.
            ipo) SPR_FLAGS+=" -ipo -ffat-lto-objects"; SPR_LDFLAGS="-fuse-ld=lld" ;;
            *) echo "ERROR: unknown SPR_TUNE=$SPR_TUNE" >&2; exit 2 ;;
        esac
        ;;
    *) echo "ERROR: unknown SPR_COMPILER=$SPR_COMPILER" >&2; exit 2 ;;
esac
module load $MPI_MODULE
module load metis/5.1.0 fftw3/3.3.10 boost/1.84.0

export CC CXX I_MPI_CC=$CC I_MPI_CXX=$CXX OMPI_CC=$CC OMPI_CXX=$CXX
export FOAM_EXTRA_CFLAGS="$SPR_FLAGS"
export FOAM_EXTRA_CXXFLAGS="$SPR_FLAGS"

# Command-line settings override etc/prefs.sh; the shared tree's prefs are not modified.
source "$OPENFOAM_DIR/etc/bashrc" \
    WM_COMPILER_TYPE=system WM_COMPILER=$WM_COMPILER WM_CC=$CC WM_CXX=$CXX \
    WM_MPLIB=$WM_MPLIB_SEL WM_COMPILE_OPTION=OptSPR$SPR_TUNE$OPT_SUFFIX

# OpenFOAM's environment loader is not errexit/nounset-safe; go strict afterwards.
set -eEo pipefail
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR

# Read-only site dependency tree. Its PT-Scotch is OpenMPI-linked: keep it out of Intel MPI builds.
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512
if [[ $SPR_MPI == intelmpi ]]; then
    export PTSCOTCH_ARCH_PATH=none
else
    # Gadi's mpicc gives no -L path; point OpenFOAM at the Open MPI prefix explicitly.
    mpi_prefix=$(ompi_info --path prefix --parsable); export MPI_ARCH_PATH=${mpi_prefix#path:prefix:} FOAM_MPI=sys-openmpi
fi
export WM_NCOMPPROCS=$BUILD_JOBS
export I_MPI_LIBRARY_KIND=release
export FOAM_EXTRA_CFLAGS="$SPR_FLAGS"
export FOAM_EXTRA_CXXFLAGS="$SPR_FLAGS"
export FOAM_EXTRA_LDFLAGS="${SPR_LDFLAGS:-}"
# Object/dependency files stay on the build node, leaving scratch inodes free.
export FOAM_BUILDROOT=${PBS_JOBFS:?Submit through PBS with jobfs}/$SPR_COMPILER$SPR_TUNE$OPT_SUFFIX-objects
mkdir "$FOAM_BUILDROOT"
export FOAM_MODULE_PREFIX=none

[[ $WM_PROJECT_VERSION == v2512 && $WM_PRECISION_OPTION == DP && $WM_LABEL_SIZE == 32 ]]
[[ $WM_OPTIONS == linux64${WM_COMPILER}DPInt32OptSPR$SPR_TUNE$OPT_SUFFIX ]]
[[ $WM_MPLIB == "$WM_MPLIB_SEL" ]]
[[ $(git -C "$OPENFOAM_DIR" rev-parse HEAD) == "$EXPECTED_COMMIT" ]]
[[ ! -e $WM_PROJECT_DIR/platforms/$WM_OPTIONS ]] || {
    echo "ERROR: $WM_PROJECT_DIR/platforms/$WM_OPTIONS already exists" >&2
    exit 1
}

. "$WM_PROJECT_DIR/wmake/scripts/have_scotch"
have_scotch
if [[ $SPR_MPI == intelmpi ]] && have_ptscotch; then
    echo "ERROR: OpenMPI PT-Scotch must not enter the Intel MPI build" >&2
    exit 2
fi

{
    echo "source_commit=$(git -C "$OPENFOAM_DIR" rev-parse HEAD)"
    echo "spr_compiler=$SPR_COMPILER"
    echo "spr_tune=${SPR_TUNE:-base}"
    echo "spr_mpi=$SPR_MPI ($MPI_MODULE)"
    echo "wm_options=$WM_OPTIONS"
    echo "wm_mplib=$WM_MPLIB"
    echo "foam_mpi=$FOAM_MPI"
    echo "compiler=$($CXX --version | head -1)"
    echo "assembler=$($($CC -print-prog-name=as) --version 2>/dev/null | head -1)"
    echo "mpi=$(mpiexec --version 2>&1 | head -1)"
    echo "foam_extra_cxxflags=$FOAM_EXTRA_CXXFLAGS"
    echo "foam_extra_ldflags=$FOAM_EXTRA_LDFLAGS"
    echo "cxxflags=$(wmake -show-cxxflags)"
    echo "build_jobs=$BUILD_JOBS"
    echo "foam_buildroot=$FOAM_BUILDROOT"
    echo "host_cpu=$(grep -m1 'model name' /proc/cpuinfo)"
    echo "timestamp_utc=$(date -u +%FT%TZ)"
} | tee "$LOG_DIR/build-info-${PBS_JOBID:-manual}.txt"
module list > "$LOG_DIR/modules-${PBS_JOBID:-manual}.txt" 2>&1

cd "$WM_PROJECT_DIR"
./Allwmake -s -j "$BUILD_JOBS" -q

for executable in simpleFoam potentialFoam decomposePar renumberMesh foamDictionary; do
    executable_path=$(command -v "$executable")
    [[ $executable_path == "$WM_PROJECT_DIR/platforms/$WM_OPTIONS/bin/$executable" ]]
    "$executable" -help >/dev/null 2>&1
    echo "$executable=$executable_path"
done

ldd "$(command -v simpleFoam)" > "$LOG_DIR/simpleFoam-ldd-${PBS_JOBID:-manual}.txt"
! grep -q 'not found' "$LOG_DIR/simpleFoam-ldd-${PBS_JOBID:-manual}.txt"
readelf -p .comment "$(command -v simpleFoam)" > "$LOG_DIR/simpleFoam-comment-${PBS_JOBID:-manual}.txt"
echo "BUILD_SUCCESS compiler=$SPR_COMPILER tune=${SPR_TUNE:-base} mpi=$SPR_MPI options=$WM_OPTIONS"
