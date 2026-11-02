#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=4,mem=16gb,walltime=00:20:00,storage=scratch/ya41
#PBS -N gcc_identity
#PBS -j oe

# Box identity test of the Bonus components built with GCC 14.2.0 -O3 -march=sapphirerapids (build_tools_gcc_jobfs.sh)
# against stock simpleFoam of the GCC OpenFOAM build (linux64GccDPInt32OptSPR). GCC fuses multiply-adds across
# statements by default (-ffp-contract=fast), so bit-identity of the replaced kernels must be shown for this compiler.
# Skewed 3D box, serial and 4 ranks, shipped fvSolution.fixedIter, fvSchemes, controlDict.noWrite (endTime 20) and
# function objects. Arms:
#   basic         stock simpleFoam (GCC build)
#   bonus_v3      preload compactAddr:fixedIterSolvers:fusedVectorSolve:stepReuse (GCC), STEP_REUSE_OFF=wss,fco, simpleFoamFused
#   bonus_v3_off  bonus_v3 with STEP_REUSE_OFF=wss,fco,fixed,gcache,pcg (hooks call the stock functions)
# usage: qsub -v WORK=<new dir>,TIMERDC=<dir with the GCC coded timer: build_gcc_OptSPR_<job>/dynamicCode_gcc> identity_test_gcc.sh
set -euo pipefail
: "${WORK:?qsub -v WORK=<dir>}"
O=/scratch/ya41/kn3896/openfoam
T=$O/tools
SHIP=/scratch/ya41/kn3896/openfoam/contest-results/compliance/reference/run/case/occDrivAerStaticMesh/system
TIMER=${TIMERDC:?qsub -v TIMERDC=...}
G=linux64GccDPInt32OptSPR
CALIB=$T/compactAddr/lib/$G/libcompactAddr.so
FIS=$T/fixedIterSolvers/lib/$G/libfixedIterSolvers.so
FLIB=$T/fusedVectorSolve/lib/$G/libfusedVectorSolve.so
RLIB=$T/stepReuse/lib/$G/libstepReuse.so
SFF=$T/simpleFoamFused/bin/$G/simpleFoamFused
set +eu
module purge; module load gcc/14.2.0 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0 > /dev/null 2>&1
source $O/OpenFOAM-v2512/etc/bashrc WM_COMPILER_TYPE=system WM_COMPILER=Gcc \
    WM_CC=gcc WM_CXX=g++ WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR > /dev/null 2>&1
set -eu
[[ $WM_OPTIONS == "$G" ]]
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 I_MPI_FABRICS=shm FOAM_SIGFPE=true \
    STEP_REUSE_REPORT=1 FUSED_VECTOR_SOLVE_REPORT=1 SIMPLEFOAM_FUSED_REPORT=1 BOX_SKEW=1
mkdir -p "$WORK" && cd "$WORK"
for l in $CALIB $FIS $FLIB $RLIB $SFF $(command -v simpleFoam); do echo "sha256 $(sha256sum "$l" | cut -d' ' -f1) $l"; done
fail=()

ARMS="basic bonus_v3 bonus_v3_off"
for a in $ARMS; do
    C=$WORK/$a
    bash "$T/fusedVectorSolve/test/make_box_case.sh" "$C" > /dev/null
    cp "$SHIP/fvSolution.fixedIter" "$C/system/fvSolution"
    cp "$SHIP/fvSchemes" "$C/system/fvSchemes"
    sed -E 's/^endTime(\s+)[0-9]+;/endTime\120;/' "$SHIP/controlDict.noWrite" > "$C/system/controlDict"
    cp -r "$SHIP/probes" "$C/system/"
    for f in FOwallClockTimeStatistics fieldMinMax solverInfo forceCoeffsAll calc_cp wallShearStress; do cp "$SHIP/$f" "$C/system/"; done
    mkdir -p "$C/system/include"
    sed -E 's/^(fixedWallPatches|forcesWallPatches)(\s+).*/\1\2ground;/' "$SHIP/include/caseDefinition" > "$C/system/include/caseDefinition"
    mkdir -p "$C/dynamicCode/platforms/$WM_OPTIONS"
    cp -r "$TIMER/wallClockTimeStatistics" "$C/dynamicCode/"
    cp -r "$TIMER/platforms/$WM_OPTIONS/lib" "$C/dynamicCode/platforms/$WM_OPTIONS/"
    pre="" off="" app=simpleFoam
    case $a in
    bonus_v3)     pre=$CALIB:$FIS:$FLIB:$RLIB off=wss,fco app=$SFF ;;
    bonus_v3_off) pre=$CALIB:$FIS:$FLIB:$RLIB off=wss,fco,fixed,gcache,pcg app=$SFF ;;
    esac
    if true; then
        cmp "$SHIP/fvSolution.fixedIter" "$C/system/fvSolution"
        [[ $(diff "$SHIP/controlDict.noWrite" "$C/system/controlDict" | grep -c '^[<>]') == 2 ]]
    fi
    (cd "$C" && blockMesh > log.blockMesh 2>&1 && cp -r . "${C}_par")
    (cd "$C" && LD_PRELOAD=$pre STEP_REUSE_OFF=$off $app > log.serial 2>&1) || fail+=("$a serial exit $?")
    (cd "${C}_par" && decomposePar -force > log.decompose 2>&1 \
        && mpirun -np 4 -genv LD_PRELOAD "$pre" -genv STEP_REUSE_OFF "$off" $app -parallel > log.par 2>&1) || fail+=("$a par exit $?")
    for m in serial par; do
        d=$C; [[ $m == par ]] && d=${C}_par
        echo "== $a $m: exit-ok steps=$(grep -c '^Time =' "$d/log.$m") End=$(grep -cx End "$d/log.$m")"
        grep -hE '^(fixedIterSolvers|stepReuse|compactAddr|fusedVectorSolve|simpleFoamFused):' "$d/log.$m" | sort -u | sed 's/^/    /' || true
    done
done

# --- comparisons against basic ---------------------------------------------------------
types() {   # normalised line types of a log (case dir -> CASE, numbers -> #), minus the application line
    python3 - "$1" "$2" <<'PY'
import re, sys
num = re.compile(r'[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?')
seen = set()
for l in open(sys.argv[1], errors='replace'):
    l = ' '.join(num.sub('#', l.replace(sys.argv[2], 'CASE')).split())
    if l and not l.startswith('Exec ') and l not in seen:
        seen.add(l); print(l)
PY
}
for a in bonus_v3 bonus_v3_off; do
    for m in serial par; do
        s=$WORK/basic f=$WORK/$a; [[ $m == par ]] && s=${s}_par f=${f}_par
        r=""
        diff -q <(grep -E 'Solving for|continuity errors' "$s/log.$m") <(grep -E 'Solving for|continuity errors' "$f/log.$m") > /dev/null || r="$r residuals"
        n=0
        while IFS= read -r p; do
            n=$((n + 1))
            cmp -s "$s/$p" "$f/$p" || r="$r $p"
        done < <(cd "$s" && find postProcessing -type f -not -path '*/wallClockTimeStatistics/*' | sort)
        [[ $(grep -c '^Time =' "$s/log.$m") == 20 && $(grep -c '^Time =' "$f/log.$m") == 20 ]] || r="$r steps"
        missing=$(comm -23 <(types "$s/log.$m" "$s" | sort) <(types "$f/log.$m" "$f" | sort))
        echo "IDENTITY $a $m: ${r:-residuals and all $n postProcessing files bit-identical to basic}"
        if [[ -n $missing ]]; then echo "F5 $a $m: basic line types missing:"; sed 's/^/    /' <<< "$missing"; else echo "F5 $a $m: every basic line type present"; fi
        [[ -z $r ]] || fail+=("$a $m identity:$r")
        [[ -z $missing ]] || fail+=("$a $m F5")
    done
done
# --- report lines (master) -----------------------------------------------------------------
need() { grep -qF "$2" "$WORK/$1/log.serial" && grep -qF "$2" "$WORK/${1}_par/log.par" || fail+=("$1 missing report: $2"); }
never() { ! grep -qF "$2" "$WORK/$1/log.serial" "$WORK/${1}_par/log.par" || fail+=("$1 unexpected report: $2"); }
for f in k omega; do
    need bonus_v3 "fixedIterSolvers: smoothSolver $f: fixed-iteration path"
    need bonus_v3_off "fixedIterSolvers: smoothSolver $f: stock solve"
done
for g in 'grad(p)' 'grad(k)' 'grad(omega)'; do need bonus_v3 "fixedIterSolvers: cache: $g cached in code"; done
need bonus_v3 'fixedIterSolvers: cache: grad(p) stored when asked again'
# (with pcgFixedPath.C in the build add: need bonus_v3 'residual norm after the last iteration only';
#  need bonus_v3_off 'fixedIterSolvers: PCG '; never bonus_v3_off 'fixed iteration count')
need bonus_v3 'fusedVectorSolve: U fused path (v3 folded copy)'
need bonus_v3 'compactAddr: GaussSeidel: 16-bit neighbour offsets, folded source copy'
never bonus_v3 'stock solve'
never bonus_v3_off 'cached in code'
never bonus_v3_off 'fixed-iteration path'
for a in bonus_v3 bonus_v3_off; do
    never $a 'stepReuse: wallShearStress: evaluated in write'
    never $a 'stepReuse: forceCoeffs: coefficient fields only with writeFields'
done

if [[ ${#fail[@]} == 0 ]]; then echo "IDENTITY_ALL_PASS"; else printf 'IDENTITY_FAIL: %s\n' "${fail[@]}"; exit 1; fi
