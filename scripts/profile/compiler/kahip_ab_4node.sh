#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:25:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N kahip_ab_4n
#PBS -j oe

# KaHIP vs scotch at 4 nodes with the current Bonus setup + UCX_ZCOPY_THRESH=inf UCX_RNDV_THRESH=inf.
#   scotch - prepared case decomp-screen/case416_w05 (team scotch (4 104)/(1 0.5), renumbered)
#   kahip  - team (an0642) KaHIP multiLevel partition, fields regenerated (kahip/regen_bench_181335722/case_kahip;
#            Oct job 181335722: -0.36% with the setup of then). Processor dirs are only linked, never written.
# 1) one gamgCensus run on kahip (10 steps): per-rank compute/neighbours and the addressing statistics (shows
#    whether the KaHIP processor meshes were renumbered); scotch census: job 182165954 (base arm)
# 2) A/B, 20 steps, 3 repeats, alternating start (scotch kahip | kahip scotch | scotch kahip)
# Different partitions are not bit-identical: Cd at step 20 is recorded for each run.

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
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/kahip/ab_4node_${PBS_JOBID%%.*}
KAH=/scratch/ya41/kn3896/openfoam/contest-results/kahip/regen_bench_181335722/case_kahip
CLIB=/scratch/ya41/kn3896/openfoam/tools/gamgCensus/lib/linux64IcxDPInt32OptSPR/libgamgCensus.so
REFLOG=/scratch/ya41/kn3896/openfoam/contest-results/diag/vtune_ranks_4node_182139223/simpleFoam.log
STEPS=20
CASEDIR=$ROOT/case
CASEK=$ROOT/case_kahip
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "KAHIP_AB_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CUR/logFiles/$l.log"; done
[[ -f $LIB && -f $FLIB && -f $RLIB && -f $CLIB && -x $SFF && -f $REFLOG ]]
for l in 30_potentialFoam 40_applyBoundaryLayer; do grep -q "^End" "$KAH/logFiles/$l.log"; done
[[ $(ls -d "$KAH"/processors416_* | wc -l) == 4 ]]
for l in "$LIB" "$FLIB" "$RLIB" "$CLIB" "$SFF"; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
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
# kahip: same system/constant/mesh, KaHIP processor dirs
mkdir -p "$CASEK/constant"
cp -a "$CASEDIR/system" "$CASEDIR/dynamicCode" "$CASEK/"
cp -a "$CUR"/constant/{transportProperties,turbulenceProperties} "$CASEK/constant/"
ln -s "$MESH" "$CASEK/constant/polyMesh"
for p in "$KAH"/processors416_*; do ln -s "$p" "$CASEK/${p##*/}"; done

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

for c in "$CASEDIR" "$CASEK"; do
    cp "$TEAM/dicts/4node/controlDict" "$c/system/controlDict"
    sed -i "s/^endTime .*/endTime           $STEPS;/" "$c/system/controlDict"
    echo "libs (\"$LIB\");" >> "$c/system/controlDict"
done

pres() { awk -v n="$2" '/^Time = /{t++} t>n{exit} /Solving for p,/{print}' "$1"; }
ucxT=(-genv UCX_ZCOPY_THRESH inf -genv UCX_RNDV_THRESH inf)
declare -A CASES=([scotch]="$CASEDIR" [kahip]="$CASEK")
run() {   # $1 arm (scotch|kahip), $2 repeat or "census"
    local a=$1 r=$2 dir="$ROOT/$1_r$2" c=${CASES[$1]} steps=$STEPS pre=(-genv LD_PRELOAD "$FLIB:$RLIB")
    mkdir -p "$dir"
    cd "$c"
    if [[ $r == census ]]; then
        steps=10; mkdir -p "$dir/census"
        pre=(-genv LD_PRELOAD "$CLIB:$FLIB:$RLIB" -genv STEP_REUSE_CENSUS "$dir/census")
    fi
    sed -i "s/^endTime .*/endTime           $steps;/" system/controlDict
    rm -rf postProcessing
    echo "RUN_START arm=$a repeat=$r utc=$(date -u +%FT%TZ)"
    local status=0
    I_MPI_PIN_PROCESSOR_LIST=$PIN_NORMAL timeout --kill-after=30s 600s \
        mpirun -np 416 -ppn 104 -hosts "$HOSTS_NORMAL" "${ucxT[@]}" "${glibcT[@]}" "${pre[@]}" \
        -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1 -genv SIMPLEFOAM_FUSED_REPORT 1 \
        "$SFF" -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    local cd=NA
    if [[ -d postProcessing ]]; then
        f=$(find postProcessing -name coefficient.dat | head -1)
        [[ -n $f ]] && cd=$(awk '!/^#/{l=$2} END{print l}' "$f")
        n=$(find postProcessing | wc -l); tar czf "$dir/postProcessing.tar.gz" postProcessing
        [[ $(tar tzf "$dir/postProcessing.tar.gz" | wc -l) == "$n" ]] && rm -rf postProcessing
    fi
    [[ -z $(ls -d "$c"/processors416_*/[1-9]* 2>/dev/null) ]]   # nothing written
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $steps ]]
    grep -qF 'simpleFoamFused: pressure equation steps fused' "$dir/simpleFoam.log"
    grep -qF 'stepReuse: laplacian: fused Gauss linear assembly' "$dir/simpleFoam.log"
    local same=NA
    if [[ $a == scotch && $r != census ]]; then
        same=DIFFER; cmp -s <(pres "$dir/simpleFoam.log" 10) <(pres "$REFLOG" 10) && same=IDENTICAL
    fi
    echo "RUN_SUCCESS arm=$a repeat=$r Cd_step$steps=$cd p_residuals_vs_ref=$same $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
    cd "$ROOT"
}

ucx_info -v | head -1
run kahip census
rep=0
for order in "scotch kahip" "kahip scotch" "scotch kahip"; do
    rep=$((rep + 1))
    for a in $order; do run "$a" "$rep"; done
done

python3 - "$ROOT" <<'PY'
import re, sys, statistics as st
from pathlib import Path
root = Path(sys.argv[1])
t = {}
for f in sorted(root.glob('*_r[0-9]*/simpleFoam.log')):
    a, r = f.parent.name.rsplit('_r', 1)
    v = [float(l.split('=')[1]) for l in f.read_text().splitlines() if 'Average wall-clock time per time step' in l]
    t.setdefault(a, []).append(v[-1])
b = st.mean(t['scotch'])
print('| Arm | repeats (s/step) | mean | min | vs base (mean) | vs base (min) |')
print('|---|---|---:|---:|---:|---:|')
for a in ['scotch', 'kahip']:
    v = t[a]
    print(f'| {a} | {" / ".join(f"{x:.4f}" for x in v)} | {st.mean(v):.4f} | {min(v):.4f} | '
          f'{100 * (st.mean(v) / b - 1):+.2f}% | {100 * (min(v) / min(t["scotch"]) - 1):+.2f}% |')
PY
# pack the linked case copy (inode quota)
for d in case case_kahip; do
    [[ -d $d ]] || continue
    n=$(find "$d" | wc -l); tar czf "$d.tar.gz" -C "$(dirname "$d")" "$(basename "$d")"
    [[ $(tar tzf "$d.tar.gz" | wc -l) == "$n" ]] && rm -rf "$d"
done
echo "KAHIP_AB_COMPLETE root=$ROOT"
