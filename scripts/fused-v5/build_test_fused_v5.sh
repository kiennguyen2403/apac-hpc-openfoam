#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=16,mem=64gb
#PBS -l walltime=00:40:00
#PBS -l jobfs=60GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N fused_v5
#PBS -j oe

# fusedVectorSolve v5 (2026-11-03): the fused U solve now computes every per-iteration convergence check (residual of
# each component with its global sum after every nSweeps sweeps, as stock smoothSolver::solve) and runs at log 2,
# printing those checks in stock's format. This job
#   1. builds libfusedVectorSolve for ICX (production, -xSAPPHIRERAPIDS -fp-model=precise) and GCC (-march=sapphirerapids)
#      from a jobfs copy of our tree (src/ and applications/ are archived); the v3 builds are in lib_v3/;
#   2. log-2 proof on the 3D box (tools/fusedVectorSolve/test/identity_test.sh, serial and 4 ranks), U with its own
#      solver dict plus 'log 2;', stock simpleFoam vs the preloaded library: fields, Solving-for lines, Normalisation
#      factor lines and every 'smoothSolver:  Iteration N residual = r' line must match:
#        shipped   the shipped U controls (maxIter 15, nSweeps 1)
#        nsweeps2  nSweeps 2 (checks at 2, 4, ..., 14; 16 sweeps)
#        miniter   minIter 20 > maxIter 15 (no Iteration 0 line, no checks printed past maxIter, 20 sweeps)
#        stack     shipped controls, full Bonus preload + simpleFoamFused
#   3. full-stack identity at log 1 with the shipped dicts and function objects (tools/fixedIterSolvers/test/identity_test.sh).

set -eEo pipefail
: "${PBS_JOBID:?submit through PBS}"
O=/scratch/ya41/kn3896/openfoam
K=$O/OpenFOAM-v2512
S=$O/tools/fusedVectorSolve
OUT=$O/contest-results/compliance/fused_v5_${PBS_JOBID%%.*}
JT=$PBS_JOBFS/OpenFOAM-v2512
mkdir -p "$OUT"
exec > >(tee "$OUT/job.log") 2>&1
cp "$0" "$OUT/build_test_fused_v5.sh" 2>/dev/null || true
cp "$S/fusedSolveSegregated.C" "$OUT/fusedSolveSegregated.C.v5"

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
    rm -f "$S/lib/$OPT/libfusedVectorSolve.so"
    rm -rf "$S/Make/$OPT"
    ( cd "$S" && wmake libso ) > "$OUT/wmake_$1.log" 2>&1 || { tail -30 "$OUT/wmake_$1.log"; exit 1; }
    [[ -f $S/lib/$OPT/libfusedVectorSolve.so ]]
    grep -q 'per-iteration residual check (v5)' "$S/lib/$OPT/libfusedVectorSolve.so"
    echo "BUILT $1 $(sha256sum "$S/lib/$OPT/libfusedVectorSolve.so" | cut -c1-16) $S/lib/$OPT/libfusedVectorSolve.so"
}
( build icx )
( build gcc )

# --- 2. log-2 proof --------------------------------------------------------------------------------------------------
I=linux64IcxDPInt32OptSPR
FLIB=$S/lib/$I/libfusedVectorSolve.so
STACK=$O/tools/compactAddr/lib/$I/libcompactAddr.so:$O/tools/fixedIterSolvers/lib/$I/libfixedIterSolvers.so:$FLIB:$O/tools/stepReuse/lib/$I/libstepReuse.so
SFF=$O/tools/simpleFoamFused/bin/$I/simpleFoamFused
BASE='solver smoothSolver; smoother GaussSeidel; tolerance 0; relTol 0;'
fail=()
log2() {   # $1 name, $2 U solver entries, $3 preload, [$4 app]
    local W=$PBS_JOBFS/log2_$1
    echo "=== log2 $1: U { $2 }"
    WORK=$W U_DICT="$2" PRELOAD=$3 APP=${4:-} bash "$S/test/identity_test.sh" | tee "$OUT/log2_$1.txt" || { fail+=("$1: test exit"); return; }
    grep -q 'per-iteration residual check (v5)' "$W/fused/log.serial" || fail+=("$1: fused path not used")
    [[ $(grep -c 'IDENTITY .*: fields all identical; solver lines identical' "$OUT/log2_$1.txt") == 2 ]] || fail+=("$1: mismatch")
    mkdir -p "$OUT/log2_$1"
    for d in stock fused stock_par fused_par; do
        m=serial; [[ $d == *_par ]] && m=par
        cp "$W/$d/log.$m" "$OUT/log2_$1/log.$d"
    done
    # first step's U lines of both arms, for the report
    awk '/^Time = 1$/ {p=1} /^Time = 2$/ {exit} p && /Normalisation|Iteration|Solving for U/' "$W/stock/log.serial" > "$OUT/log2_$1/step1_U_stock.txt"
    awk '/^Time = 1$/ {p=1} /^Time = 2$/ {exit} p && /Normalisation|Iteration|Solving for U/' "$W/fused/log.serial" > "$OUT/log2_$1/step1_U_fused.txt"
}
log2 shipped  "$BASE maxIter 15; nSweeps 1; log 2;"            "$FLIB"
log2 nsweeps2 "$BASE maxIter 15; nSweeps 2; log 2;"            "$FLIB"
log2 miniter  "$BASE minIter 20; maxIter 15; nSweeps 1; log 2;" "$FLIB"
log2 stack    "$BASE maxIter 15; nSweeps 1; log 2;"            "$STACK" "$SFF"

# --- 3. full-stack identity, shipped dicts, log 1 ------------------------------------------------------------------
echo "=== full-stack identity (shipped fvSolution.fixedIter, fvSchemes, controlDict.noWrite, function objects)"
WORK=$PBS_JOBFS/stack_identity bash "$O/tools/fixedIterSolvers/test/identity_test.sh" | tee "$OUT/stack_identity.txt" || fail+=("stack identity exit")
grep -q '^IDENTITY_ALL_PASS' "$OUT/stack_identity.txt" || fail+=("stack identity")

if (( ${#fail[@]} )); then printf 'FAIL: %s\n' "${fail[@]}"; echo "FUSED_V5_FAILED"; else echo "FUSED_V5_ALL_PASS"; fi
