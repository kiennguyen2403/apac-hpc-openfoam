#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=32,mem=64gb
#PBS -l walltime=01:30:00
#PBS -l jobfs=120GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N build_gcc_arch
#PBS -j oe

# OpenFOAM v2512 with GCC 14.2.0 + Intel MPI 2021.17.2, -O3 -march=$ARCH, platform linux64GccDPInt32$OPTNAME.
# usage: qsub -v ARCH=sapphirerapids,OPTNAME=OptSPR build_gcc_arch_jobfs.sh   (or ARCH=skylake-avx512,OPTNAME=OptSKX)
# FOAM_EXTRA_*FLAGS are exported AFTER sourcing etc/bashrc (bashrc unsets them; jobs 182942641/182947329 lost their
# -march that way and were compiled at GCC's default -march=broadwell), and the build stops unless wmake's compile
# line contains -march=$ARCH. Built on jobfs (our src/ and applications/ are archived; inode quota); the finished
# platform and the coded timer for it ($OUT/dynamicCode_gcc) are copied to scratch.

set -eEo pipefail
: "${PBS_JOBID:?submit through PBS}"
: "${ARCH:?qsub -v ARCH=...,OPTNAME=...}" "${OPTNAME:?qsub -v ARCH=...,OPTNAME=...}"
K=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
OUT=/scratch/ya41/kn3896/openfoam/contest-results/compiler-ab/build_gcc_${OPTNAME}_${PBS_JOBID%%.*}
LOGD=/scratch/ya41/kn3896/openfoam/logs/build/v2512-gcc-$ARCH-intelmpi
JT=$PBS_JOBFS/OpenFOAM-v2512
OPT=linux64GccDPInt32$OPTNAME
mkdir -p "$OUT" "$LOGD"
exec > >(tee "$LOGD/build-${PBS_JOBID%%.*}.log") 2>&1
trap 'echo "GCC_BUILD_FAILED line=$LINENO status=$?"' ERR
[[ ! -e $K/platforms/$OPT ]] || { echo "ERROR: $K/platforms/$OPT exists"; exit 1; }
echo "START utc=$(date -u +%FT%TZ)"

# --- tree on jobfs ------------------------------------------------------------------------
rsync -a --exclude platforms --exclude build --exclude '*.tar' --exclude tutorials --exclude 'log.*' "$K/" "$JT/"
tar -xf "$K/src.tar" -C "$JT"
tar -xf "$K/applications.tar" -C "$JT"
[[ -d $JT/src/OpenFOAM && -d $JT/applications/solvers ]]
echo "tree on jobfs: $(find "$JT" | wc -l) entries"

set +eu
module purge
module load gcc/14.2.0 binutils/2.44 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
export CC=gcc CXX=g++ I_MPI_CC=gcc I_MPI_CXX=g++

source "$JT/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Gcc WM_CC=gcc WM_CXX=g++ WM_MPLIB=INTELMPI \
    WM_COMPILE_OPTION=$OPTNAME
set -eEo pipefail
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
export FOAM_EXTRA_CFLAGS="-march=$ARCH" FOAM_EXTRA_CXXFLAGS="-march=$ARCH"   # after bashrc, which unsets them
export FOAM_BUILDROOT=$PBS_JOBFS/objects WM_NCOMPPROCS=${PBS_NCPUS:-32} I_MPI_LIBRARY_KIND=release FOAM_MODULE_PREFIX=none
mkdir -p "$FOAM_BUILDROOT"
[[ $WM_OPTIONS == "$OPT" && $WM_MPLIB == INTELMPI && $WM_PROJECT_DIR == "$JT" ]]
wmake -show-cxxflags | grep -q -- "-march=$ARCH" || { echo "ERROR: -march=$ARCH not in the compile line: $(wmake -show-cxxflags)"; exit 1; }
{
    echo "source_tree=$K (src.tar, applications.tar unpacked on jobfs)"
    echo "source_commit=$(git -C "$K" rev-parse HEAD 2>/dev/null)"
    echo "spr_compiler=gcc"; echo "spr_mpi=intelmpi (intel-mpi/2021.17.2)"; echo "wm_options=$WM_OPTIONS"
    echo "compiler=$($CXX --version | head -1)"; echo "assembler=$(as --version | head -1)"
    echo "mpi=$(mpiexec --version 2>&1 | head -1)"; echo "foam_extra_cxxflags=$FOAM_EXTRA_CXXFLAGS"
    echo "cxxflags=$(wmake -show-cxxflags)"; echo "timestamp_utc=$(date -u +%FT%TZ)"
} | tee "$LOGD/build-info-${PBS_JOBID%%.*}.txt" "$OUT/build-info.txt"

cd "$WM_PROJECT_DIR"
./Allwmake -s -j "$WM_NCOMPPROCS" -q
for e in simpleFoam potentialFoam decomposePar foamDictionary; do
    [[ $(command -v $e) == "$JT/platforms/$OPT/bin/$e" ]]; $e -help > /dev/null 2>&1
done
ldd "$(command -v simpleFoam)" > "$OUT/simpleFoam-ldd.txt"; ! grep -q 'not found' "$OUT/simpleFoam-ldd.txt"
echo "BUILT utc=$(date -u +%FT%TZ)"

# --- coded timer for this platform -------------------------------------------------------------
cp -a /scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05/dynamicCode "$OUT/dynamicCode_gcc"
( cd "$OUT/dynamicCode_gcc/wallClockTimeStatistics" && wmake libso ) > "$OUT/timer-build.log" 2>&1
ls "$OUT/dynamicCode_gcc/platforms/$OPT/lib/"libwallClockTimeStatistics_*.so
rm -rf "$OUT/dynamicCode_gcc/wallClockTimeStatistics/Make/$OPT"   # objects of the timer build (inodes)

# --- platform back into our tree -----------------------------------------------------------------
n=$(find "$JT/platforms/$OPT" | wc -l)
iu=$(lquota 2>/dev/null | awk '$1=="ya41" && $2=="scratch" {print $9}')
echo "platform entries: $n; ya41 inode usage now: ${iu:-?}"
[[ -z $iu || $(( iu + n )) -lt 197000 ]] || { echo "ERROR: inode quota too close, platform left on jobfs"; exit 1; }
rsync -a "$JT/platforms/$OPT" "$K/platforms/"
rsync -a "$JT/platforms/tools/" "$K/platforms/tools/" 2>/dev/null || true
sha256sum "$K/platforms/$OPT/bin/simpleFoam" "$K/platforms/$OPT/lib/libfiniteVolume.so" | tee -a "$OUT/build-info.txt"
echo "GCC_BUILD_SUCCESS options=$OPT utc=$(date -u +%FT%TZ)"
