#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=16,mem=64gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N rebuild_bonus_tools
#PBS -j oe

# Rebuilds the five Bonus components in the bonus folder exactly as the custom-code section of its build/build.pbs
# does (same environment: its etc/bashrc with the prefs.sh that build.pbs wrote, -O3 -xCORE-AVX512 -fp-model=precise;
# same REPO path, same sed of Make/files, same wmake calls), from the clone src/apac-hpc-openfoam after it was moved to
# the commit pinned in preprocess.sh. The OpenFOAM part of build.pbs (Allwmake) is not repeated: that platform is
# unchanged (job 182941860).

set -o pipefail
: "${PBS_JOBID:?submit through PBS}"
BON=/scratch/ya41/apac-2026-bonus/openfoam
OUT=/scratch/ya41/kn3896/openfoam/contest-results/bonus-folder/rebuild_tools_${PBS_JOBID%%.*}
mkdir -p "$OUT"
exec > >(tee "$OUT/build.log") 2>&1

module purge
module load intel-mpi/2021.17.2
module load intel-compiler-llvm/2025.3.2
cd "$BON/build"
source "../src/OpenFOAM-v2512/etc/bashrc"
[[ $WM_OPTIONS == linux64IcxDPInt32Opt ]] || { echo "ERROR WM_OPTIONS=$WM_OPTIONS"; exit 1; }
wmake -show-cxxflags | tee "$OUT/cxxflags.txt" | grep -q -- -xCORE-AVX512 || true
grep -q -- -xCORE-AVX512 "$OUT/cxxflags.txt" || { echo "ERROR: prefs.sh flags missing"; exit 1; }
echo "clone HEAD (logical path): $(git -C "$BON/src/apac-hpc-openfoam" rev-parse --short HEAD)"
echo "pinned in preprocess.sh:   $(grep -oE 'checkout [0-9a-f]+' "$BON/build/preprocess.sh")"

# --- verbatim from build/build.pbs ----------------------------------------------------------------
REPO=$WM_PROJECT_DIR/../apac-hpc-openfoam
LIBDIR=$REPO/platforms/$WM_OPTIONS/lib
BINDIR=$REPO/platforms/$WM_OPTIONS/bin
mkdir -p $LIBDIR $BINDIR

cd $REPO/tools/compactAddr
sed -i "s#^LIB = .*#LIB = $LIBDIR/libcompactAddr#" Make/files
wmake libso

cd $REPO/tools/fixedIterSolvers
sed -i "s#^LIB = .*#LIB = $LIBDIR/libfixedIterSolvers#" Make/files
wmake libso

cd $REPO/tools/fusedVectorSolve
sed -i "s#^LIB = .*#LIB = $LIBDIR/libfusedVectorSolve#" Make/files
wmake libso

cd $REPO/tools/stepReuse
sed -i "s#^LIB = .*#LIB = $LIBDIR/libstepReuse#" Make/files
wmake libso

cd $REPO/tools/simpleFoamFused
sed -i "s#^EXE = .*#EXE = $BINDIR/simpleFoamFused#" Make/files
wmake
# --- end of the build.pbs section ------------------------------------------------------------------

ok=1
for f in $LIBDIR/libcompactAddr.so $LIBDIR/libfixedIterSolvers.so $LIBDIR/libfusedVectorSolve.so $LIBDIR/libstepReuse.so $BINDIR/simpleFoamFused; do
    if [[ -e $f && $f -nt $OUT/cxxflags.txt ]]; then echo "BUILT $(sha256sum "$f" | cut -c1-16) $(readlink -f "$f")"; else echo "NOT REBUILT $f"; ok=0; fi
done
nm -DC $LIBDIR/libstepReuse.so > "$OUT/nm_stepReuse.txt"
if grep -qE 'wallShearStress::(execute|write)|forceCoeffs::(reset|calcForceCoeffs)' "$OUT/nm_stepReuse.txt"; then
    echo "ERROR: libstepReuse still defines a function-object replacement"; ok=0
fi
(( ok )) && echo "BONUS_TOOLS_REBUILD_OK" || echo "BONUS_TOOLS_REBUILD_FAILED"
