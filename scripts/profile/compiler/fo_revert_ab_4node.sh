#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:15:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N fo_ab_4n
#PBS -j oe

# Validation of the forward prediction in contest-results/models/REPORT.md (recorded 2026-10-25T12:00:22Z):
# cost of re-enabling the two function objects the Bonus arm defers. 4 nodes, current Bonus setup
# (compactAddr + fusedVectorSolve + stepReuse, simpleFoamFused, UCX zero-copy off), 20 steps, ONE repeat
# (user request), arm order fixed:
#   cur    - production Bonus (wallShearStress evaluated at write only, forceCoeffs fields skipped)
#   wss    - same + STEP_REUSE_OFF=wss      (stock wallShearStress every step)      predicted +0.73 .. +1.40%
#   wssfco - same + STEP_REUSE_OFF=wss,fco  (+ stock forceCoeffs fields)            predicted +1.15 .. +2.04%
# controlDict: the shipped controlDict.noWrite with endTime 20, plus the fixedIterSolvers libs entry.
# p residuals of every run (steps 1-10) must equal job 182139223's.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
RLIB=/scratch/ya41/kn3896/openfoam/tools/stepReuse/lib/linux64IcxDPInt32OptSPR/libstepReuse.so
SFF=/scratch/ya41/kn3896/openfoam/tools/simpleFoamFused/bin/linux64IcxDPInt32OptSPR/simpleFoamFused
CUR=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/models/fo_revert_ab_4node_${PBS_JOBID%%.*}
CALIB=/scratch/ya41/kn3896/openfoam/tools/compactAddr/lib/linux64IcxDPInt32OptSPR/libcompactAddr.so
CLIB=/scratch/ya41/kn3896/openfoam/tools/gamgCensus/lib/linux64IcxDPInt32OptSPR/libgamgCensus.so
REFLOG=/scratch/ya41/kn3896/openfoam/contest-results/diag/vtune_ranks_4node_182139223/simpleFoam.log
STEPS=20
CASEDIR=$ROOT/case
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "FO_AB_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CUR/logFiles/$l.log"; done
[[ -f $LIB && -f $FLIB && -f $RLIB && -f $CALIB && -x $SFF && -f $REFLOG ]]
for l in "$LIB" "$FLIB" "$RLIB" "$CALIB" "$SFF"; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
cp "$CUR/system/fvSolution.fixedIter" fvSolution.basic
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.basic > fvSolution.bonus
[[ $(grep -c 'smoothSolverFixed;' fvSolution.bonus) == 3 && $(grep -c 'PCGFixed;' fvSolution.bonus) == 1 ]]
awk '{print} /^[[:space:]]*grad\(U\);/ && incache {print "    grad(p);"; print "    grad(k);"; print "    grad(omega);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.bonus > fvSolution.tmp && mv fvSolution.tmp fvSolution.bonus
[[ $(grep -c 'grad(p);' fvSolution.bonus) == 1 && $(grep -c 'grad(omega);' fvSolution.bonus) == 1 ]]
lscpu -p=CPU,CORE,SOCKET,NODE > lscpu.txt

mkdir -p "$CASEDIR/constant"
cp -a "$CUR/system" "$CASEDIR/"
cp fvSolution.bonus "$CASEDIR/system/fvSolution"
cp -a "$CUR"/constant/{transportProperties,turbulenceProperties} "$CASEDIR/constant/"
ln -s "$MESH" "$CASEDIR/constant/polyMesh"
for p in "$CUR"/processors416_*; do ln -s "$p" "$CASEDIR/${p##*/}"; done
cp -a "$CUR/dynamicCode" "$CASEDIR/dynamicCode"

set +eu
module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
    WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
set -eu
[[ $WM_OPTIONS == linux64IcxDPInt32OptSPR ]]
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FOAM_SIGFPE=true
export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx I_MPI_DEBUG=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_EXCLUDE_LIST I_MPI_PIN_DOMAIN_SIZE
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
PIN_NORMAL=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
PIN_PERM=$(tr ',' '\n' <<< "$PIN_NORMAL" | tac | paste -sd,)
[[ $(tr ',' '\n' <<< "$PIN_NORMAL" | wc -l) == 104 && $(tr ',' '\n' <<< "$PIN_PERM" | wc -l) == 104 ]]
HOSTS_NORMAL=$(sort -u "$PBS_NODEFILE" | paste -sd,)
HOSTS_PERM=$(sort -u "$PBS_NODEFILE" | tac | paste -sd,)
echo "hosts normal: $HOSTS_NORMAL"; echo "hosts perm:   $HOSTS_PERM"
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-4n-104ppn-scatter.dat
glibcT=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)

cd "$CASEDIR"
sed -E "s/^(endTime\s+)[0-9]+;/\1$STEPS;/" system/controlDict.noWrite > system/controlDict
[[ $(diff system/controlDict.noWrite system/controlDict | grep -c '^[<>]') == 2 ]]
echo "libs (\"$LIB\");" >> system/controlDict

pres() { awk -v n="$2" '/^Time = /{t++} t>n{exit} /Solving for p,/{print}' "$1"; }
base=(-genv UCX_ZCOPY_THRESH inf -genv UCX_RNDV_THRESH inf)
declare -A OFF=(
    [cur]=""
    [wss]="wss"
    [wssfco]="wss,fco"
)
run() {   # $1 arm, $2 repeat; failures are reported, not fatal
    local a=$1 r=$2 dir="$ROOT/$1_r$2"
    mkdir -p "$dir"
    cd "$CASEDIR"
    rm -rf postProcessing
    echo "RUN_START arm=$a repeat=$r utc=$(date -u +%FT%TZ)"
    local status=0 off=()
    [[ -n ${OFF[$a]} ]] && off=(-genv STEP_REUSE_OFF "${OFF[$a]}")
    # shellcheck disable=SC2086
    I_MPI_PIN_PROCESSOR_LIST=$PIN_NORMAL timeout --kill-after=30s 600s \
        mpirun -np 416 -ppn 104 -hosts "$HOSTS_NORMAL" "${base[@]}" "${glibcT[@]}" -genv LD_PRELOAD "$CALIB:$FLIB:$RLIB" "${off[@]}" \
        -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1 -genv SIMPLEFOAM_FUSED_REPORT 1 \
        "$SFF" -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then
        n=$(find postProcessing | wc -l); tar czf "$dir/postProcessing.tar.gz" postProcessing
        [[ $(tar tzf "$dir/postProcessing.tar.gz" | wc -l) == "$n" ]] && rm -rf postProcessing
    fi
    cd "$ROOT"
    [[ -z $(ls -d "$CASEDIR"/processors416_*/[1-9]* 2>/dev/null) ]]   # nothing written (fatal if violated)
    if [[ $status != 0 ]] || ! grep -qx End "$dir/simpleFoam.log" || [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") != $STEPS ]]; then
        echo "RUN_FAILED arm=$a repeat=$r status=$status"; tail -5 "$dir/simpleFoam.log"; return 0
    fi
    local same=DIFFER
    cmp -s <(pres "$dir/simpleFoam.log" 10) <(pres "$REFLOG" 10) && same=IDENTICAL
    echo "RUN_SUCCESS arm=$a repeat=$r p_residuals=$same $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

ucx_info -v | head -1
for rep in 1; do
    for a in cur wss wssfco; do run "$a" "$rep"; done
done
for a in cur wss wssfco; do
    echo "deferral report lines, arm $a:"
    grep -hE "stepReuse: (wallShearStress|forceCoeffs)" "$ROOT/${a}_r1/simpleFoam.log" | sort | uniq -c || echo "  (none)"
done

python3 - "$ROOT" <<'PY'
import re, sys, statistics as st
from pathlib import Path
root = Path(sys.argv[1])
t = {}
for f in sorted(root.glob('*_r*/simpleFoam.log')):
    a, r = f.parent.name.rsplit('_r', 1)
    v = [float(l.split('=')[1]) for l in f.read_text().splitlines() if 'Average wall-clock time per time step' in l]
    if v:
        t.setdefault(a, []).append(v[-1])
b = st.mean(t['cur'])
pred = {'wss': (0.73, 1.40), 'wssfco': (1.15, 2.04)}
print('| Arm | s/step | vs cur | predicted | within prediction (+-0.35% noise)? |')
print('|---|---:|---:|---:|---|')
for a in ['cur', 'wss', 'wssfco']:
    v = t.get(a, [])
    if not v:
        print(f'| {a} | failed | | | |'); continue
    d = 100 * (st.mean(v) / b - 1)
    lo, hi = pred.get(a, (0, 0))
    ok = '' if a == 'cur' else ('yes' if lo - 0.35 <= d <= hi + 0.35 else 'NO')
    print(f'| {a} | {st.mean(v):.4f} | {d:+.2f}% | ' + ('-' if a == 'cur' else f'+{lo:.2f} .. +{hi:.2f}%') + f' | {ok} |')
if 'wss' in t and 'wssfco' in t:
    d = 100 * (st.mean(t['wssfco']) - st.mean(t['wss'])) / b
    print(f'\nforceCoeffs fields alone (wssfco - wss): {d:+.2f}% (predicted +0.42 .. +0.64%)')
PY
# pack the linked case copy (inode quota)
for d in case; do
    [[ -d $d ]] || continue
    n=$(find "$d" | wc -l); tar czf "$d.tar.gz" -C "$(dirname "$d")" "$(basename "$d")"
    [[ $(tar tzf "$d.tar.gz" | wc -l) == "$n" ]] && rm -rf "$d"
done
echo "FO_AB_COMPLETE root=$ROOT"
