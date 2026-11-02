#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=16,mem=64gb
#PBS -l walltime=01:00:00
#PBS -l jobfs=60GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N build_tools_gcc
#PBS -j oe

# The five Bonus components (compactAddr, fixedIterSolvers, fusedVectorSolve, stepReuse, simpleFoamFused; production
# sources = Make/files as for the compliant finals) built with GCC 14.2.0 -O3 -march=sapphirerapids against the GCC
# OpenFOAM build linux64GccDPInt32OptSPR (job 182952873). Their Make/files put the output under lib/$(WM_OPTIONS) and
# bin/$(WM_OPTIONS), so the ICX builds (lib/linux64IcxDPInt32OptSPR) are not touched.
# Our tree's src/ and applications/ are archived (inode quota), so the headers come from a copy of the tree on jobfs with
# src.tar/applications.tar unpacked and the GCC platform copied in; FOAM_EXTRA_*FLAGS are exported after etc/bashrc
# (which unsets them) and checked in the compile line.

set -eEo pipefail
: "${PBS_JOBID:?submit through PBS}"
K=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
T=/scratch/ya41/kn3896/openfoam/tools
OPT=linux64GccDPInt32OptSPR
OUT=/scratch/ya41/kn3896/openfoam/contest-results/compiler-ab/build_tools_gcc_${PBS_JOBID%%.*}
JT=$PBS_JOBFS/OpenFOAM-v2512
mkdir -p "$OUT"
exec > >(tee "$OUT/build.log") 2>&1
[[ -x $K/platforms/$OPT/bin/simpleFoam ]]

rsync -a --exclude platforms --exclude build --exclude '*.tar' --exclude tutorials --exclude 'log.*' "$K/" "$JT/"
mkdir -p "$JT/platforms"
rsync -a "$K/platforms/$OPT" "$JT/platforms/"
rsync -a "$K/platforms/tools" "$JT/platforms/"   # wmake helpers (wmkdepend, lemon) for linux64Gcc
[[ -x $JT/platforms/tools/linux64Gcc/wmkdepend ]]
tar -xf "$K/src.tar" -C "$JT"
tar -xf "$K/applications.tar" -C "$JT"

set +eu
module purge
module load gcc/14.2.0 binutils/2.44 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
export CC=gcc CXX=g++ I_MPI_CC=gcc I_MPI_CXX=g++
source "$JT/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Gcc WM_CC=gcc WM_CXX=g++ WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
set -eEo pipefail
trap 'echo "TOOLS_GCC_FAILED line=$LINENO status=$?"' ERR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none I_MPI_LIBRARY_KIND=release
export FOAM_EXTRA_CFLAGS="-march=sapphirerapids" FOAM_EXTRA_CXXFLAGS="-march=sapphirerapids"   # after bashrc
[[ $WM_OPTIONS == "$OPT" && $WM_PROJECT_DIR == "$JT" ]]
wmake -show-cxxflags | grep -q -- '-march=sapphirerapids'
echo "compiler=$(g++ --version | head -1)"
echo "cxxflags=$(wmake -show-cxxflags)"

for t in compactAddr fixedIterSolvers fusedVectorSolve stepReuse; do
    echo "== $t"
    rm -f "$T/$t/lib/$OPT/lib$t.so"
    ( cd "$T/$t" && wmake libso ) > "$OUT/$t.wmake.log" 2>&1 || { echo "wmake failed for $t:"; tail -40 "$OUT/$t.wmake.log"; exit 1; }
    grep -E ' error|Error ' "$OUT/$t.wmake.log" | head -5 || true
    [[ -f $T/$t/lib/$OPT/lib$t.so ]] || { echo "MISSING $T/$t/lib/$OPT/lib$t.so"; exit 1; }
done
echo "== simpleFoamFused"
rm -f "$T/simpleFoamFused/bin/$OPT/simpleFoamFused"
( cd "$T/simpleFoamFused" && wmake ) > "$OUT/simpleFoamFused.wmake.log" 2>&1 || { echo "wmake failed for simpleFoamFused:"; tail -40 "$OUT/simpleFoamFused.wmake.log"; exit 1; }
[[ -x $T/simpleFoamFused/bin/$OPT/simpleFoamFused ]]

# replaced symbols present (written to files: grep -q on a pipe under pipefail misreports)
nm -DC "$T/fixedIterSolvers/lib/$OPT/libfixedIterSolvers.so" > "$OUT/nm_fis.txt"
nm -DC "$T/compactAddr/lib/$OPT/libcompactAddr.so" > "$OUT/nm_calib.txt"
grep -q ' T Foam::smoothSolver::solve' "$OUT/nm_fis.txt"
grep -q ' T Foam::solution::cache' "$OUT/nm_fis.txt"
grep -q ' T Foam::lduMatrix::Amul' "$OUT/nm_calib.txt"
{
    for f in "$T"/{compactAddr,fixedIterSolvers,fusedVectorSolve,stepReuse}/lib/$OPT/*.so "$T/simpleFoamFused/bin/$OPT/simpleFoamFused"; do
        echo "$(sha256sum "$f" | cut -c1-16) $f"
    done
} | tee "$OUT/sha256.txt"
echo "TOOLS_GCC_BUILD_SUCCESS utc=$(date -u +%FT%TZ)"
