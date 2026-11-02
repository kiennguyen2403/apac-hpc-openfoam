#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=4,mem=16gb,walltime=00:20:00,storage=scratch/ya41
#PBS -N bonus_identity
#PBS -j oe

# Box identity test of the BONUS FOLDER build (/scratch/ya41/apac-2026-bonus/openfoam: build/build.pbs, prefs.sh
# -O3 -xCORE-AVX512 -fp-model=precise, platform linux64IcxDPInt32Opt in the Basic folder's tree; custom code from the
# clone src/apac-hpc-openfoam at bf45d4c). Derived from tools/fixedIterSolvers/test/identity_test.sh, which tested the
# same code under the SPR flags. Skewed 3D box, serial and 4 ranks, shipped fvSolution.fixedIter, fvSchemes,
# controlDict.noWrite (endTime 20) and function objects. Arms:
#   basic         stock simpleFoam of the bonus-folder build
#   bonus_v3      preload compactAddr:fixedIterSolvers:fusedVectorSolve:stepReuse, STEP_REUSE_OFF=wss,fco, simpleFoamFused
#   bonus_v3_off  bonus_v3 with STEP_REUSE_OFF=wss,fco,fixed,gcache,pcg (hooks call the stock functions)
# The coded timer is compiled by OpenFOAM in each case (no prebuilt copy for this platform).
# usage: qsub -v WORK=<new dir> identity_bonus_folder.sh
set -euo pipefail
: "${WORK:?qsub -v WORK=<dir>}"
O=/scratch/ya41/kn3896/openfoam
T=$O/tools
SHIP=/scratch/ya41/kn3896/openfoam/contest-results/compliance/reference/run/case/occDrivAerStaticMesh/system
TIMER=$O/contest-results/vtune/drivaer_1node_icx_180728663/case/dynamicCode
BON=/scratch/ya41/apac-2026-bonus/openfoam
set +eu
module purge; module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0 > /dev/null 2>&1
source "$BON/src/OpenFOAM-v2512/etc/bashrc" > /dev/null 2>&1
set -eu
[[ $WM_OPTIONS == linux64IcxDPInt32Opt ]]
REPO=$WM_PROJECT_DIR/../apac-hpc-openfoam
PL=$REPO/platforms/$WM_OPTIONS
CALIB=$PL/lib/libcompactAddr.so
FIS=$PL/lib/libfixedIterSolvers.so
FLIB=$PL/lib/libfusedVectorSolve.so
RLIB=$PL/lib/libstepReuse.so
SFF=$PL/bin/simpleFoamFused
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
