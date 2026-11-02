#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=32,mem=64gb
#PBS -l walltime=01:40:00
#PBS -l jobfs=60GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N build_bonus_folder
#PBS -j oe

# Runs the bonus folder's own build script, /scratch/ya41/apac-2026-bonus/openfoam/build/build.pbs, unchanged
# (its prefs.sh: Icx, Intel MPI, -O3 -xCORE-AVX512 -fp-model=precise, WM_COMPILE_OPTION Opt; OpenFOAM built in
# ../src/OpenFOAM-v2512, a symlink to the Basic folder's tree; then compactAddr, fixedIterSolvers, fusedVectorSolve,
# stepReuse and simpleFoamFused from the clone ../src/apac-hpc-openfoam at bf45d4c). The only difference from
# "qsub build.pbs": FOAM_BUILDROOT on the node's jobfs, so the ~20k object/dependency files do not land on scratch
# (ya41 inode quota). Object location does not change the binaries. Then checks what was built.

set -o pipefail
: "${PBS_JOBID:?submit through PBS}"
BD=/scratch/ya41/apac-2026-bonus/openfoam/build
LOG=/scratch/ya41/kn3896/openfoam/logs/build/bonus-folder
mkdir -p "$LOG"
exec > >(tee "$LOG/build-${PBS_JOBID%%.*}.log") 2>&1
echo "BUILD_START utc=$(date -u +%FT%TZ) build.pbs sha256=$(sha256sum "$BD/build.pbs" | cut -d' ' -f1)"

export FOAM_BUILDROOT=$PBS_JOBFS/objects
mkdir -p "$FOAM_BUILDROOT"
cd "$BD"
bash ./build.pbs
echo "BUILD_PBS_EXIT=$? utc=$(date -u +%FT%TZ)"

# --- checks -------------------------------------------------------------------------
set +eu
module purge > /dev/null 2>&1
module load intel-mpi/2021.17.2 intel-compiler-llvm/2025.3.2 > /dev/null 2>&1
source "$BD/../src/OpenFOAM-v2512/etc/bashrc" > /dev/null 2>&1
REPO=$WM_PROJECT_DIR/../apac-hpc-openfoam
P=$REPO/platforms/$WM_OPTIONS
ok=1
echo "WM_PROJECT_DIR=$WM_PROJECT_DIR (resolved $(readlink -f "$WM_PROJECT_DIR")) WM_OPTIONS=$WM_OPTIONS"
echo "cxxflags=$(wmake -show-cxxflags 2>/dev/null)"
echo "clone HEAD=$(git -C "$REPO" rev-parse --short HEAD) dirty files: $(git -C "$REPO" status --short | wc -l)"
for f in "$WM_PROJECT_DIR/platforms/$WM_OPTIONS/bin/simpleFoam" "$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib/libfiniteVolume.so" \
         "$P/lib/libcompactAddr.so" "$P/lib/libfixedIterSolvers.so" "$P/lib/libfusedVectorSolve.so" "$P/lib/libstepReuse.so" \
         "$P/bin/simpleFoamFused"; do
    if [[ -e $f ]]; then echo "BUILT $(sha256sum "$f" | cut -c1-16) $f"; else echo "MISSING $f"; ok=0; fi
done
for s in "smoothSolver::solve" "solution::cache"; do nm -DC "$P/lib/libfixedIterSolvers.so" 2>/dev/null | grep -q " T Foam::$s" || { echo "MISSING symbol $s"; ok=0; }; done
ldd "$P/bin/simpleFoamFused" | grep 'not found' && ok=0
echo "lnInclude dirs created: $(find "$WM_PROJECT_DIR/src" -maxdepth 4 -type d -name lnInclude | wc -l)"
lquota 2>/dev/null | grep ya41
if (( ok )); then echo "BONUS_FOLDER_BUILD_OK"; else echo "BONUS_FOLDER_BUILD_INCOMPLETE"; fi
