#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=16,mem=64gb
#PBS -l walltime=00:40:00
#PBS -l jobfs=60GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N build_stepreuse
#PBS -j oe

# Rebuild of libstepReuse after removing wallShearStressLazy.C and forceCoeffsLazy.C from its Make/files (2026-11-03),
# for both compilers it exists for:
#   icx  ICX 2025.3.2 -O3 -xSAPPHIRERAPIDS -fp-model=precise against linux64IcxDPInt32OptSPR (production)
#   gcc  GCC 14.2.0   -O3 -march=sapphirerapids             against linux64GccDPInt32OptSPR
# Our tree's src/ and applications/ are archived (inode quota), so the headers come from a jobfs copy of the tree with
# the sources unpacked and the platforms + wmake helpers copied in. FOAM_EXTRA_*FLAGS are exported after etc/bashrc
# (which unsets them) and checked in the compile line. The previous ICX build is lib_v4/libstepReuse.so.

set -eEo pipefail
: "${PBS_JOBID:?submit through PBS}"
K=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
S=/scratch/ya41/kn3896/openfoam/tools/stepReuse
OUT=/scratch/ya41/kn3896/openfoam/contest-results/compliance/stepreuse_rebuild_${PBS_JOBID%%.*}
JT=$PBS_JOBFS/OpenFOAM-v2512
mkdir -p "$OUT"
exec > >(tee "$OUT/build.log") 2>&1
! grep -qE 'wallShearStressLazy|forceCoeffsLazy' "$S/Make/files"

rsync -a --exclude platforms --exclude build --exclude '*.tar' --exclude tutorials --exclude 'log.*' "$K/" "$JT/"
mkdir -p "$JT/platforms"
rsync -a "$K/platforms/linux64IcxDPInt32OptSPR" "$K/platforms/linux64GccDPInt32OptSPR" "$K/platforms/tools" "$JT/platforms/"
tar -xf "$K/src.tar" -C "$JT"
tar -xf "$K/applications.tar" -C "$JT"

build() {   # $1 icx|gcc ; in a subshell
    set +eu
    module purge
    if [[ $1 == icx ]]; then
        module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
        export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
        source "$JT/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
        FLAGS="-xSAPPHIRERAPIDS -fp-model=precise"; OPT=linux64IcxDPInt32OptSPR
    else
        module load gcc/14.2.0 binutils/2.44 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
        export CC=gcc CXX=g++ I_MPI_CC=gcc I_MPI_CXX=g++
        source "$JT/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Gcc WM_CC=gcc WM_CXX=g++ WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
        FLAGS="-march=sapphirerapids"; OPT=linux64GccDPInt32OptSPR
    fi
    set -eEo pipefail
    export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
    export FOAM_EXTRA_CFLAGS="$FLAGS" FOAM_EXTRA_CXXFLAGS="$FLAGS"   # after bashrc
    [[ $WM_OPTIONS == "$OPT" && $WM_PROJECT_DIR == "$JT" ]]
    wmake -show-cxxflags > "$OUT/cxxflags_$1.txt"
    for f in $FLAGS; do grep -q -- "$f" "$OUT/cxxflags_$1.txt"; done
    rm -f "$S/lib/$OPT/libstepReuse.so"
    rm -rf "$S/Make/$OPT"
    ( cd "$S" && wmake libso ) > "$OUT/wmake_$1.log" 2>&1 || { tail -30 "$OUT/wmake_$1.log"; exit 1; }
    [[ -f $S/lib/$OPT/libstepReuse.so ]]
    nm -DC "$S/lib/$OPT/libstepReuse.so" > "$OUT/nm_$1.txt"
    if grep -qE 'wallShearStress::(execute|write)|forceCoeffs::(reset|calcForceCoeffs)' "$OUT/nm_$1.txt"; then
        echo "ERROR: $1 build still defines a function-object replacement"; exit 1
    fi
    echo "BUILT $1 $(sha256sum "$S/lib/$OPT/libstepReuse.so" | cut -c1-16) $S/lib/$OPT/libstepReuse.so"
}
( build icx )
( build gcc )
echo "STEPREUSE_REBUILD_SUCCESS utc=$(date -u +%FT%TZ)"
