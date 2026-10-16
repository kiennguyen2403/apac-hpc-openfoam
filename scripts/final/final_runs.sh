#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N final_runs
#PBS -j oe

# Final runs, occDrivAerStaticMesh 65M. Submit per node count, e.g.
#   qsub -l ncpus=416,mem=400gb -v NNODES=4,REPEATS=3 final_runs.sh
#   qsub -l ncpus=208,mem=400gb -v NNODES=2,REPEATS=1 final_runs.sh
#   qsub -l ncpus=104,mem=400gb -v NNODES=1,REPEATS=1 final_runs.sh
# Optional: ARMS="bonus" (default "basic bonus"); REF_CD=<Cd of a basic run> to check a bonus-only job.
# Steps: the team's controlDict for that node count, unchanged (1n 40, 2n 80, 4n 150).
# Two arms per repeat, interleaved on the same nodes (basic, bonus, basic, bonus, ...):
#   basic : Basic tier - team fvSolution (stock solvers), team decomposition (scotch; 2n/4n
#           domains (N 104), domainWeights (1 0.5)), ICX OptSPR + Intel MPI, scatter pinning,
#           tuned Allreduce file, glibcT malloc environment
#   bonus : same + tools/fixedIterSolvers (smoothSolverFixed / PCGFixed via fvSolution solver
#           names and controlDict libs) + tools/fusedVectorSolve (fused Ux/Uy/Uz GaussSeidel,
#           mpirun -genv LD_PRELOAD) + "grad(p);" in the fvSolution cache block; bit-identical results,
#           so Cd must match basic exactly.
#           Since job 181376137 (-4.31% at 4n, bit-identical) also: tools/stepReuse (LD_PRELOAD:
#           GaussSeidel folded per-sweep copy, kOmegaSST shared divU/convection, wallShearStress
#           evaluated at write only), fixedIterSolvers v2, fusedVectorSolve v3, and "grad(k); grad(omega);"
#           in the cache block.
#           (Job 181332282 ran bonus without the fused vector solve.)
# Prepared cases (prep chain done once, reused): 1n vtune/drivaer_1node_icx_180728663/case,
# 2n final/case208 (parallel prep done here on first use), 4n decomp-screen/case416_w05.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
: "${NNODES:?qsub -v NNODES=1|2|4}"
REPEATS=${REPEATS:-1}
ARMS=${ARMS:-basic bonus}
REF_CD=${REF_CD:-}
NP=$((NNODES * 104))
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
RES=/scratch/ya41/kn3896/openfoam/contest-results
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
RLIB=/scratch/ya41/kn3896/openfoam/tools/stepReuse/lib/linux64IcxDPInt32OptSPR/libstepReuse.so
TIMER_SRC=$RES/vtune/drivaer_1node_icx_180728663/case/dynamicCode
ROOT=$RES/final/runs_${NNODES}n_${PBS_JOBID%%.*}
case $NNODES in
    1) CASE=$RES/vtune/drivaer_1node_icx_180728663/case; LOGS=$CASE/../logFiles; fh=(-fileHandler hostCollated); STEPS=40 ;;
    2) CASE=$RES/final/case208; LOGS=$CASE/logFiles; fh=(-fileHandler collated -ioRanks '(0 104)'); STEPS=80 ;;
    4) CASE=$RES/decomp-screen/case416_w05; LOGS=$CASE/logFiles; fh=(-fileHandler collated -ioRanks '(0 104 208 312)'); STEPS=150 ;;
    *) echo "NNODES must be 1, 2 or 4"; exit 1 ;;
esac

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "FINAL_FAILED line=$LINENO status=$?"; [[ -f $ROOT/fvSolution.basic ]] && cp "$ROOT/fvSolution.basic" "$CASE/system/fvSolution"; cp "$TEAM/dicts/${NNODES}node/controlDict" "$CASE/system/controlDict" 2>/dev/null || true' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == $NNODES ]]
cp "$TEAM/run/${NNODES}node/intel-scatter-cpus.awk" "$TEAM/run/${NNODES}node/intelmpi-allreduce-${NNODES}n-104ppn-scatter.dat" .
[[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' "$TEAM/dicts/${NNODES}node/controlDict") == $STEPS ]]
[[ -f $LIB && -f $FLIB && -f $RLIB ]]
for a in $ARMS; do [[ $a == basic || $a == bonus ]]; done
grep -q '^End' "$LOGS/10_decomposePar.log"

# --- parallel prep (2-node case on first use; team chain, site module) -------------------
prep_done=1
for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$LOGS/$l.log" 2>/dev/null || prep_done=0; done
if [[ $prep_done == 0 ]]; then
(
    set +eu
    module purge
    module load OpenFOAM/2512
    set -eu
    cd "$CASE"
    par=(mpirun -np $NP --map-by ppr:104:node --bind-to core)
    echo "PREP_START utc=$(date -u +%FT%TZ)"
    "${par[@]}" renumberMesh -constant -overwrite -parallel "${fh[@]}" > logFiles/20_renumberMesh.log 2>&1
    "${par[@]}" potentialFoam -initialiseUBCs -parallel "${fh[@]}" > logFiles/30_potentialFoam.log 2>&1
    "${par[@]}" applyBoundaryLayer -ybl "0.0450244" -parallel "${fh[@]}" > logFiles/40_applyBoundaryLayer.log 2>&1
    for l in 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' logFiles/$l.log; done
    echo "PREP_DONE utc=$(date -u +%FT%TZ)"
)
fi
[[ -d $CASE/dynamicCode ]] || cp -a "$TIMER_SRC" "$CASE/dynamicCode"
ls "$CASE/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
cp "$CASE/system/fvSolution" fvSolution.basic
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.basic > fvSolution.bonus
[[ $(grep -c 'smoothSolverFixed;' fvSolution.bonus) == 3 && $(grep -c 'PCGFixed;' fvSolution.bonus) == 1 ]]
# Bonus also caches grad(p) (job 181366105: -2.30%, bit-identical): reused once per step while p is unchanged
awk '{print} /^[[:space:]]*grad\(U\);/ && incache {print "    grad(p);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.bonus > fvSolution.bonus.tmp && mv fvSolution.bonus.tmp fvSolution.bonus
[[ $(grep -c 'grad(p);' fvSolution.bonus) == 1 ]]
# ... and grad(k), grad(omega): computed for CDkOmega, reused by the corrected laplacian terms
awk '{print} /^[[:space:]]*grad\(p\);/ && incache {print "    grad(k);"; print "    grad(omega);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.bonus > fvSolution.bonus.tmp && mv fvSolution.bonus.tmp fvSolution.bonus
[[ $(grep -c 'grad(k);' fvSolution.bonus) == 1 && $(grep -c 'grad(omega);' fvSolution.bonus) == 1 ]]
if [[ -d $CASE/postProcessing ]]; then mv "$CASE/postProcessing" "$ROOT/postProcessing.before"; fi

# --- solver environment ------------------------------------------------------------
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
I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
export I_MPI_PIN_PROCESSOR_LIST
[[ $(tr ',' '\n' <<< "$I_MPI_PIN_PROCESSOR_LIST" | wc -l) == 104 ]]
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-${NNODES}n-104ppn-scatter.dat
glibcT=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)
{
    echo "simpleFoam=$(command -v simpleFoam)"
    echo "compiler=$(icpx --version | head -1)"
    echo "mpi=$(mpirun --version 2>&1 | head -1)"
    echo "plugin=$LIB sha256=$(sha256sum "$LIB" | cut -d' ' -f1)"
    echo "fused=$FLIB sha256=$(sha256sum "$FLIB" | cut -d' ' -f1)"
    echo "stepReuse=$RLIB sha256=$(sha256sum "$RLIB" | cut -d' ' -f1)"
} > build-metadata.txt
mpirun -np $NP -ppn 104 hostname | sort | uniq -c > placement.txt
if [[ $(wc -l < placement.txt) != $NNODES ]] || awk '$1 != 104 {bad=1} END {exit !bad}' placement.txt; then
    echo "PLACEMENT_WRONG:"; cat placement.txt; exit 1
fi

run() {   # $1 arm (basic|bonus), $2 repeat
    local a=$1 r=$2 dir=$ROOT/run_$1_$2 status=0
    mkdir "$dir"
    cd "$CASE"
    cp "$TEAM/dicts/${NNODES}node/controlDict" system/controlDict
    cp "$ROOT/fvSolution.$a" system/fvSolution
    local pre=()
    if [[ $a == bonus ]]; then
        echo "libs (\"$LIB\");" >> system/controlDict
        pre=(-genv LD_PRELOAD "$FLIB:$RLIB" -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1)
    fi
    cp system/controlDict system/fvSolution "$dir/"
    echo "RUN_START arm=$a repeat=$r utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 2400s mpirun -np $NP -ppn 104 "${glibcT[@]}" "${pre[@]}" simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cp "$ROOT/fvSolution.basic" system/fvSolution
    cp "$TEAM/dicts/${NNODES}node/controlDict" system/controlDict
    cd "$ROOT"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    if [[ $a == bonus ]]; then
        grep -q 'fusedVectorSolve: U fused path (v3 folded copy)' "$dir/simpleFoam.log"
        grep -q 'stepReuse: GaussSeidel: folded' "$dir/simpleFoam.log"
        grep -q 'stepReuse: kOmegaSST: shared' "$dir/simpleFoam.log"
        grep -q 'stepReuse: wallShearStress: evaluated in write' "$dir/simpleFoam.log"
    fi
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS arm=$a repeat=$r Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
for ((r = 1; r <= REPEATS; r++)); do for a in $ARMS; do run $a $r; done; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"   # case restored

python3 - "$ROOT" "$NNODES" "$STEPS" "$REPEATS" "$ARMS" "$REF_CD" <<'PY'
import re, sys
from pathlib import Path
root, nn, steps, reps = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
arms, ref_cd = sys.argv[5].split(), (sys.argv[6] if len(sys.argv) > 6 else '')
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
def t(a, r): return float(pat.findall((root / f'run_{a}_{r}' / 'simpleFoam.log').read_text())[-1])
def cd(a, r):
    rows = [l.split() for p in (root / f'run_{a}_{r}' / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return rows[-1][1]
both = 'basic' in arms and 'bonus' in arms
out = [f'# Final runs, {nn} node(s) x 104 ranks, occDrivAerStaticMesh 65M, {steps} steps, {reps} repeat(s), arms: {" ".join(arms)}', '',
       'basic = Basic tier; bonus = basic + fixedIterSolvers plugin + fused vector solve + stepReuse (LD_PRELOAD) + grad(p), grad(k), grad(omega) cache. Arms interleaved on the same nodes.', '',
       '| Repeat | ' + ' | '.join(f'{a} s/step' for a in arms) + (' | bonus vs basic' if both else '') + ' | ' + ' | '.join(f'Cd {a} (step {steps})' for a in arms) + ' | Cd check |',
       '|---:|' + '---:|' * (len(arms) * 2 + (1 if both else 0)) + ':---|']
for r in range(1, reps + 1):
    ts = [t(a, r) for a in arms]; cds = [cd(a, r) for a in arms]
    ref = cd('basic', r) if both else ref_cd
    check = ('identical' if all(c == ref for c in cds) else 'DIFFERENT') + ('' if both else f' vs REF_CD {ref}' if ref else '')
    out.append(f'| {r} | ' + ' | '.join(f'{x:.4f}' for x in ts) + (f' | {100 * (t("bonus", r) / t("basic", r) - 1):+.2f}%' if both else '')
               + ' | ' + ' | '.join(cds) + f' | {check if ref else "n/a"} |')
best = {a: min(t(a, r) for r in range(1, reps + 1)) for a in arms}
out += [''] + [f'Best {a}: **{v:.4f} s/step**' for a, v in best.items()]
(root / 'REPORT.md').write_text('\n'.join(out) + '\n'); print('\n'.join(out))
PY
echo "FINAL_COMPLETE root=$ROOT"
