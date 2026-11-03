#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=4,mem=16gb,walltime=00:15:00,storage=scratch/ya41
#PBS -N fused_identity
#PBS -j oe

# Bit-identity test of libfusedVectorSolve against stock OpenFOAM on the 3D box case
# (make_box_case.sh), serial and 4 ranks. usage: qsub -v WORK=<empty scratch dir>[,PRELOAD=<lib[:lib]>] identity_test.sh
# PRELOAD defaults to libfusedVectorSolve; any preloaded library must leave results bit-identical.
# Optional for the test arm: SOLVER_SED (sed expression applied to system/fvSolution, e.g.
# 's/solver GAMG;/solver GAMGFused;/') and LIBS (library path(s) added to controlDict libs).
# Optional for both arms: CACHE_KOMEGA=1 (cache grad(k), grad(omega)), WSS=1 (wallShearStress
# function object on the wall, written at the last step), FCO=1 (forceCoeffs on the wall every step,
# coefficient.dat compared), DEBUG_CACHE=1 (log solution cache use),
# BOX_SKEW=1 (non-orthogonal box, read by make_box_case.sh).
# Optional for the test arm: APP=<solver executable> instead of simpleFoam (e.g. simpleFoamFused).
# Optional for both arms: U_DICT='<solver entries>' gives U its own solver dict (k and omega keep
# theirs), e.g. with 'log 2;' to print every per-iteration check; the logs' Normalisation factor and
# Iteration lines are then compared too.
set -euo pipefail
: "${WORK:?qsub -v WORK=<dir>}"
T=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve
FL=${PRELOAD:-$T/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so}
set +eu
module purge; module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 > /dev/null 2>&1
source /scratch/ya41/kn3896/openfoam/OpenFOAM-v2512/etc/bashrc WM_COMPILER_TYPE=system WM_COMPILER=Icx \
    WM_CC=icx WM_CXX=icpx WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR > /dev/null 2>&1
set -eu
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 I_MPI_FABRICS=shm FOAM_SIGFPE=true FUSED_VECTOR_SOLVE_REPORT=1 CSR_AMUL_REPORT=1 GAMG_FUSED_REPORT=1 STEP_REUSE_REPORT=1 SIMPLEFOAM_FUSED_REPORT=1
mkdir -p "$WORK" && cd "$WORK"
for l in ${FL//:/ }; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
for v in stock fused; do
    bash "$T/test/make_box_case.sh" "$WORK/$v" > /dev/null
    if [[ -n ${U_DICT:-} ]]; then sed -i "s/^    \"(U|k|omega)\" { \(.*\) }\$/    U { $U_DICT }\n    \"(k|omega)\" { \1 }/" "$WORK/$v/system/fvSolution"; grep -q "^    U { $U_DICT }" "$WORK/$v/system/fvSolution"; fi
    pre=""; [[ $v == fused ]] && pre=$FL
    app=simpleFoam; [[ $v == fused && -n ${APP:-} ]] && app=$APP
    [[ ${CACHE_KOMEGA:-0} == 1 ]] && echo 'cache { grad(k); grad(omega); }' >> "$WORK/$v/system/fvSolution"
    [[ ${DEBUG_CACHE:-0} == 1 ]] && echo 'DebugSwitches { solution 1; }' >> "$WORK/$v/system/controlDict"
    fos=""
    [[ ${WSS:-0} == 1 ]] && fos="$fos wss { type wallShearStress; libs (fieldFunctionObjects); patches (ground); writeControl writeTime; }"
    [[ ${FCO:-0} == 1 ]] && fos="$fos fco { type forceCoeffs; libs (forces); patches (ground); rho rhoInf; rhoInf 1; CofR (0 0 0); liftDir (0 0 1); dragDir (1 0 0); pitchAxis (0 1 0); magUInf 10; lRef 1; Aref 1; writeControl timeStep; writeInterval 1; }"
    if [[ -n $fos ]]; then echo "functions { $fos }" >> "$WORK/$v/system/controlDict"; fi
    if [[ $v == fused && -n ${SOLVER_SED:-} ]]; then sed -i "$SOLVER_SED" "$WORK/$v/system/fvSolution"; grep -E 'solver' "$WORK/$v/system/fvSolution" | tr -s ' ' | head -3; fi
    if [[ $v == fused && -n ${LIBS:-} ]]; then echo "libs ($(for l in ${LIBS//:/ }; do printf '"%s" ' "$l"; done));" >> "$WORK/$v/system/controlDict"; tail -1 "$WORK/$v/system/controlDict"; fi
    (cd "$v" && blockMesh > log.blockMesh 2>&1 && cp -r . "../${v}_par"
     LD_PRELOAD=$pre $app > log.serial 2>&1
     echo "$v serial exit=$? steps=$(grep -c '^Time =' log.serial) $(grep -E 'fusedVectorSolve|csrAmul|gamgFused|stepReuse|simpleFoamFused' log.serial | sort -u | tr '\n' ' ' || true)")
    (cd "${v}_par" && decomposePar -force > log.decompose 2>&1 && mpirun -np 4 -genv LD_PRELOAD "$pre" $app -parallel > log.par 2>&1
     echo "$v parallel exit=$? steps=$(grep -c '^Time =' log.par) $(grep -E 'fusedVectorSolve|csrAmul|gamgFused|stepReuse|simpleFoamFused' log.par | sort -u | tr '\n' ' ' || true)"
     reconstructPar -latestTime > log.recon 2>&1)
done
for m in serial par; do
    s=stock; f=fused; [[ $m == par ]] && s=stock_par f=fused_par
    r=""; for fld in U p k omega nut; do cmp -s "$s/20/$fld" "$f/20/$fld" || r="$r $fld"; done
    if [[ ${FCO:-0} == 1 ]]; then
        cmp -s <(grep -hv '^#' "$s"/postProcessing/fco/*/coefficient*.dat) <(grep -hv '^#' "$f"/postProcessing/fco/*/coefficient*.dat) || r="$r coefficient.dat"
        [[ -s $(ls "$s"/postProcessing/fco/*/coefficient*.dat | head -1) ]] || r="$r coefficient.dat-missing"
    fi
    if [[ ${WSS:-0} == 1 ]]; then
        cmp -s "$s/20/wallShearStress" "$f/20/wallShearStress" || r="$r wallShearStress"
        diff -q <(cat "$s"/postProcessing/wss/*/*.dat) <(cat "$f"/postProcessing/wss/*/*.dat) > /dev/null || r="$r wss.dat"
        [[ -s "$s/20/wallShearStress" ]] || r="$r wallShearStress-missing"
    fi
    pat='Solving for|continuity'
    if [[ -n ${U_DICT:-} ]]; then
        pat="$pat|Normalisation factor|:  Iteration [0-9]+ residual"
        echo "LOG2 $m: Iteration lines stock $(grep -cE ':  Iteration [0-9]+ residual' "$s/log.$m") fused $(grep -cE ':  Iteration [0-9]+ residual' "$f/log.$m"); Normalisation lines stock $(grep -c 'Normalisation factor' "$s/log.$m") fused $(grep -c 'Normalisation factor' "$f/log.$m")"
    fi
    if diff -q <(grep -E "$pat" "$s/log.$m") <(grep -E "$pat" "$f/log.$m") > /dev/null; then l=identical; else l=DIFFER; fi
    echo "IDENTITY $m: fields ${r:-all identical}; solver lines $l"
done
