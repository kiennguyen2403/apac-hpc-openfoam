#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:15:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N pcg_ab_4n
#PBS -j oe

# A/B of the coarsest-level PCG fixed-count path (fixedIterSolvers v3b, pcgFixedPath.C): the coarsest GAMG level
# runs PCG with minIter 35 = maxIter 35, so the residual norm after each of the first 34 iterations (one Allreduce
# each, 408 per step over 12 V-cycles) only feeds a convergence test that cannot change the count; the path skips it.
# Bit-identical (box identity jobs 182885471/182885540). Expected 0.2-0.9% at 4 nodes (models/REPORT.md, E1).
# 4 nodes, compliant Bonus setup (shipped fvSolution.fixedIter + controlDict.noWrite with endTime 20, no libs line,
# LD_PRELOAD compactAddr:fixedIterSolvers:fusedVectorSolve:stepReuse, simpleFoamFused), 20 steps, REPEATS (default 1),
# arm order fixed:
#   cur  - STEP_REUSE_OFF=wss,fco,pcg  (= compliant finals 182882563)
#   pcg  - STEP_REUSE_OFF=wss,fco
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
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/pcg/ab_4node_${PBS_JOBID%%.*}
REPEATS=${REPEATS:-1}
CALIB=/scratch/ya41/kn3896/openfoam/tools/compactAddr/lib/linux64IcxDPInt32OptSPR/libcompactAddr.so
CLIB=/scratch/ya41/kn3896/openfoam/tools/gamgCensus/lib/linux64IcxDPInt32OptSPR/libgamgCensus.so
REFLOG=/scratch/ya41/kn3896/openfoam/contest-results/diag/vtune_ranks_4node_182139223/simpleFoam.log
STEPS=20
CASEDIR=$ROOT/case
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "PCG_AB_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CUR/logFiles/$l.log"; done
[[ -f $LIB && -f $FLIB && -f $RLIB && -f $CALIB && -x $SFF && -f $REFLOG ]]
for l in "$LIB" "$FLIB" "$RLIB" "$CALIB" "$SFF"; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
SHIP=$TEAM/run/case/occDrivAerStaticMesh/system
cmp "$SHIP/fvSolution.fixedIter" "$CUR/system/fvSolution.fixedIter"
lscpu -p=CPU,CORE,SOCKET,NODE > lscpu.txt

mkdir -p "$CASEDIR/constant"
cp -a "$CUR/system" "$CASEDIR/"
cp "$SHIP/fvSolution.fixedIter" "$CASEDIR/system/fvSolution"
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

pres() { awk -v n="$2" '/^Time = /{t++} t>n{exit} /Solving for p,/{print}' "$1"; }
base=(-genv UCX_ZCOPY_THRESH inf -genv UCX_RNDV_THRESH inf)
declare -A OFF=(
    [cur]="wss,fco,pcg"
    [pcg]="wss,fco"
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
        mpirun -np 416 -ppn 104 -hosts "$HOSTS_NORMAL" "${base[@]}" "${glibcT[@]}" -genv LD_PRELOAD "$CALIB:$LIB:$FLIB:$RLIB" "${off[@]}" \
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
for ((rep = 1; rep <= REPEATS; rep++)); do
    for a in cur pcg; do run "$a" "$rep"; done
done
for a in cur pcg; do
    echo "fixedIterSolvers report lines, arm $a (repeat 1):"
    grep -hE "fixedIterSolvers: (PCG|smoothSolver)" "$ROOT/${a}_r1/simpleFoam.log" | sort | uniq -c || echo "  (none)"
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
print('| Arm | repeats (s/step) | mean | min | vs cur (mean) | vs cur (min) |')
print('|---|---|---:|---:|---:|---:|')
for a in ['cur', 'pcg']:
    v = t.get(a, [])
    if not v:
        print(f'| {a} | failed | | | | |'); continue
    print(f'| {a} | {" / ".join(f"{x:.4f}" for x in v)} | {st.mean(v):.4f} | {min(v):.4f} | '
          f'{100 * (st.mean(v) / b - 1):+.2f}% | {100 * (min(v) / min(t["cur"]) - 1):+.2f}% |')
print('\npredicted (models/REPORT.md E1): -0.2 .. -0.9%; paired noise about 0.35%')
PY
# pack the linked case copy (inode quota)
for d in case; do
    [[ -d $d ]] || continue
    n=$(find "$d" | wc -l); tar czf "$d.tar.gz" -C "$(dirname "$d")" "$(basename "$d")"
    [[ $(tar tzf "$d.tar.gz" | wc -l) == "$n" ]] && rm -rf "$d"
done
echo "PCG_AB_COMPLETE root=$ROOT"
