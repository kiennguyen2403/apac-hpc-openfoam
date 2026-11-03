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
# Bonus switches (one per item put to the judges; defaults = what is switched on now):
#   GCACHE=1      grad(p), grad(k), grad(omega) cached in code (libfixedIterSolvers; 0 -> STEP_REUSE_OFF+=gcache)
#   FIXED_ITER=0  default since 2026-11-03 (every convergence check must be computed; more cheaply is allowed):
#                 STEP_REUSE_OFF+=fixed, k and omega on the stock smoothSolver (the log must show it). FIXED_ITER=1
#                 (fixed-iteration smoothSolver path of libfixedIterSolvers: residual after the last sweep only)
#                 skips those checks: not allowed in a scored run (an unscored "how fast it could go" slide only).
#   FUSED_U=1     libfusedVectorSolve preloaded, which must be v5 (fused Ux/Uy/Uz solve computing every per-iteration
#                 check, log-2 proof job 182994209; the log must show its v5 line); 0 -> not preloaded (U on the
#                 stock smoothSolver). Versions up to v3 skipped the checks.
#   Function objects always run exactly as in the stock code: the wallShearStress and forceCoeffs deferrals were
#   removed from libstepReuse on 2026-11-03 (before that they were switched off with STEP_REUSE_OFF=wss,fco).
# Steps (Addendum 2026-09-25, F1): the shipped controlDict.noWrite with only endTime changed
# (1n 40, 2n 80, 4n 150); fvSolution.fixedIter, fvSchemes and controlDict.noWrite must match the
# shipped case, and caseDefinition may differ from it only in nCores / decompositionMethod /
# nHierarchical (scotch domains/domainWeights live in decomposeParDict scotchCoeffs).
# Both arms run the shipped fvSolution and the run controlDict unchanged (no libs line) since 2026-11-02.
# (Jobs before 2026-10-24 used the team's controlDict, which also set runTimeModifiable no; Bonus jobs before
# 2026-11-02 renamed solvers and added cache entries in fvSolution and loaded fixedIterSolvers v2 through controlDict
# libs, and deferred wallShearStress to write time and skipped the forceCoeffs fields.)
# Two arms per repeat, interleaved on the same nodes (basic, bonus, basic, bonus, ...):
#   basic : Basic tier - team fvSolution (stock solvers), team decomposition (scotch; 2n/4n
#           domains (N 104), domainWeights (1 0.5)), ICX OptSPR + Intel MPI, scatter pinning,
#           tuned Allreduce file, glibcT malloc environment, and (since job 182166167: -2.7% at 4n,
#           bit-identical, transport only) UCX_ZCOPY_THRESH=inf UCX_RNDV_THRESH=inf: inter-node messages
#           go through UCX's pre-registered buffers instead of registering solver buffers (zero-copy /
#           rendezvous), which slowed the compute of ranks with off-node neighbours (~5 us per off-node
#           face per step; jobs 182165344, 182165954)
#   bonus : same + tools/fixedIterSolvers (until 2026-11-02: smoothSolverFixed / PCGFixed via fvSolution
#           solver names and controlDict libs) + tools/fusedVectorSolve (fused Ux/Uy/Uz GaussSeidel,
#           mpirun -genv LD_PRELOAD) + grad(p) cached (until 2026-11-02 "grad(p);" in the fvSolution cache
#           block); bit-identical results, so Cd must match basic exactly.
#           Since job 181376137 (-4.31% at 4n, bit-identical) also: tools/stepReuse (LD_PRELOAD:
#           GaussSeidel folded per-sweep copy, kOmegaSST shared divU/convection, wallShearStress
#           evaluated at write only), fixedIterSolvers v2, fusedVectorSolve v3, and "grad(k); grad(omega);"
#           in the cache block. Since job 181379590 (-3.53% at 4n, bit-identical) stepReuse also fuses
#           divDevReff and the kOmegaSST per-cell maths (stepReuse lib_v1 = the build measured in 181376794).
#           Since job 181412786 (-8.78% at 4n, bit-identical) stepReuse also fuses the Gauss laplacian,
#           the linearUpwindV convection of U and the Gauss gradient, and skips the forceCoeffs coefficient
#           fields without writeFields (stepReuse lib_v2 = the build measured at 1.4242, job 181390148).
#           Since job 182135600 (-2.72% at 4n, bit-identical) the bonus arm runs tools/simpleFoamFused
#           (simpleFoam with fused momentum/pressure steps) instead of simpleFoam; basic runs simpleFoam.
#           Since job 182137125 (-0.37% at 4n, within noise, bit-identical) stepReuse also fuses the
#           cellLimited gradient of U (stepReuse lib_v3 = the build measured at 1.2675, job 182135721).
#           Since job 182180991 (-1.6% at 4n, bit-identical; box identity job 182180957) the bonus arm also
#           preloads tools/compactAddr first: GaussSeidel and Amul read uint16 neighbour offsets and uint8
#           faces-per-cell instead of int32 addressing (the time loop is memory-bandwidth bound, job 182171696).
#           Since 2026-11-02 (compliance rework, box identity job 182881133) fixedIterSolvers v3 is preloaded
#           instead: it replaces the stock smoothSolver::solve (fixed-iteration path when the shipped controls
#           make early convergence impossible, stock solve otherwise) and solution::cache (adds grad(p), grad(k),
#           grad(omega)), so fvSolution and controlDict stay as shipped; PCGFixed dropped (never engaged: the
#           coarsest-level PCG has a non-zero default tolerance). LD_PRELOAD order:
#           compactAddr : fixedIterSolvers : fusedVectorSolve : stepReuse; STEP_REUSE_OFF=wss,fco.
#           (A PCG::scalarSolve replacement for the coarsest GAMG level, minIter 35 = maxIter 35 so the residual
#           norm after each of the first 34 iterations is skipped, was bit-identical (box jobs 182885471/182885540)
#           but +0.02% at 4 nodes (job 182885925, noise): not built, source kept as pcgFixedPath.C.)
#           (Job 181332282 ran bonus without the fused vector solve.)
# Prepared cases (prep chain done once, reused): 1n vtune/drivaer_1node_icx_180728663/case,
# 2n final/case208 (parallel prep done here on first use), 4n decomp-screen/case416_w05.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
: "${NNODES:?qsub -v NNODES=1|2|4}"
REPEATS=${REPEATS:-1}
ARMS=${ARMS:-basic bonus}
REF_CD=${REF_CD:-}
GCACHE=${GCACHE:-1}
FIXED_ITER=${FIXED_ITER:-0}
FUSED_U=${FUSED_U:-1}
NP=$((NNODES * 104))
TEAM=/scratch/ya41/kn3896/openfoam/contest-results/compliance/reference   # read-only copy of the team files (the team folder moved)
SHIP=$TEAM/run/case/occDrivAerStaticMesh
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
RES=/scratch/ya41/kn3896/openfoam/contest-results
# COMPILER=icx (default; ICX 2025.3.2 -O3 -xSAPPHIRERAPIDS -fp-model=precise, platform linux64IcxDPInt32OptSPR) or
# COMPILER=gcc (GCC 14.2.0 -O3 -march=sapphirerapids, platform linux64GccDPInt32OptSPR, build job 182952873; Bonus
# components built with GCC by build_tools_gcc_jobfs.sh). Both arms use the chosen compiler's build.
COMPILER=${COMPILER:-icx}
case $COMPILER in
    icx) OPT=linux64IcxDPInt32OptSPR ;;
    gcc) OPT=linux64GccDPInt32OptSPR
         GCCTIMER=${GCCTIMER:-/scratch/ya41/kn3896/openfoam/contest-results/compiler-ab/build_gcc_OptSPR_182952873/dynamicCode_gcc} ;;
    *) echo "COMPILER must be icx or gcc"; exit 1 ;;
esac
FIS=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/$OPT/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/$OPT/libfusedVectorSolve.so
RLIB=/scratch/ya41/kn3896/openfoam/tools/stepReuse/lib/$OPT/libstepReuse.so
SFF=/scratch/ya41/kn3896/openfoam/tools/simpleFoamFused/bin/$OPT/simpleFoamFused
CALIB=/scratch/ya41/kn3896/openfoam/tools/compactAddr/lib/$OPT/libcompactAddr.so
TIMER_SRC=$RES/vtune/drivaer_1node_icx_180728663/case/dynamicCode
PFX=""; if [[ $COMPILER == gcc ]]; then PFX=gcc_; fi
ROOT=$RES/final/runs_${PFX}${NNODES}n_${PBS_JOBID%%.*}
case $NNODES in
    1) CASE=$RES/vtune/drivaer_1node_icx_180728663/case; LOGS=$CASE/../logFiles; fh=(-fileHandler hostCollated); STEPS=40 ;;
    2) CASE=$RES/final/case208; LOGS=$CASE/logFiles; fh=(-fileHandler collated -ioRanks '(0 104)'); STEPS=80 ;;
    4) CASE=$RES/decomp-screen/case416_w05; LOGS=$CASE/logFiles; fh=(-fileHandler collated -ioRanks '(0 104 208 312)'); STEPS=150 ;;
    *) echo "NNODES must be 1, 2 or 4"; exit 1 ;;
esac

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "FINAL_FAILED line=$LINENO status=$?"; cp "$ROOT/controlDict.run" "$CASE/system/controlDict" 2>/dev/null || true' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == $NNODES ]]
cp "$TEAM/run/${NNODES}node/intel-scatter-cpus.awk" "$TEAM/run/${NNODES}node/intelmpi-allreduce-${NNODES}n-104ppn-scatter.dat" .
# Case files against the shipped case (Addendum F1)
for f in controlDict.noWrite fvSolution.fixedIter fvSchemes; do cmp "$SHIP/system/$f" "$CASE/system/$f"; done
bad=$(diff "$SHIP/system/include/caseDefinition" "$CASE/system/include/caseDefinition" | grep '^[<>]' | grep -vE '^[<>] (nCores|decompositionMethod|nHierarchical)\s' || true)
[[ -z $bad ]] || { echo "caseDefinition differs from the shipped case outside nCores/decompositionMethod/nHierarchical:"; echo "$bad"; false; }
sed -E "s/^(endTime\s+)[0-9]+;/\1$STEPS;/" "$SHIP/system/controlDict.noWrite" > controlDict.run
[[ $(diff "$SHIP/system/controlDict.noWrite" controlDict.run | grep -c '^[<>]') == 2 ]]
for f in include/caseDefinition decomposeParDict; do echo "== $f: shipped vs case"; diff "$SHIP/system/$f" "$CASE/system/$f" || true; done > compliance.txt
echo "== controlDict: shipped controlDict.noWrite vs run" >> compliance.txt
diff "$SHIP/system/controlDict.noWrite" controlDict.run >> compliance.txt || true
[[ -f $FIS && -f $FLIB && -f $RLIB && -f $CALIB && -x $SFF ]]
for a in $ARMS; do [[ $a == basic || $a == bonus ]]; done
for s in GCACHE FIXED_ITER FUSED_U; do [[ ${!s} == 0 || ${!s} == 1 ]]; done
# Bonus preload list and switched-off parts
BONUS_OFF=none   # switches nothing off (no STEP_REUSE_OFF key is a substring of "none")
[[ $GCACHE == 1 ]] || BONUS_OFF+=,gcache
[[ $FIXED_ITER == 1 ]] || BONUS_OFF+=,fixed
if [[ $FUSED_U == 1 ]]; then BONUS_PRELOAD=$CALIB:$FIS:$FLIB:$RLIB; else BONUS_PRELOAD=$CALIB:$FIS:$RLIB; fi
{
    echo "== bonus arm: $SFF -parallel, shipped fvSolution and run controlDict unchanged (no libs line)"
    echo "   GCACHE=$GCACHE FIXED_ITER=$FIXED_ITER FUSED_U=$FUSED_U"
    echo "   LD_PRELOAD=$BONUS_PRELOAD"
    echo "   STEP_REUSE_OFF=$BONUS_OFF"
} >> compliance.txt
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
if [[ $COMPILER == gcc && ! -d $CASE/dynamicCode/platforms/$OPT/lib ]]; then
    # the coded timer compiled for the GCC platform (same source, SHA1 6e8d888...; no source tree here to compile it)
    grep -q 6e8d888ad28997c1bc9eebd1c08a1a15774516a4 "$CASE/dynamicCode/wallClockTimeStatistics/Make/files"
    mkdir -p "$CASE/dynamicCode/platforms/$OPT"
    cp -a "$GCCTIMER/platforms/$OPT/lib" "$CASE/dynamicCode/platforms/$OPT/"
fi
ls "$CASE/dynamicCode/platforms/$OPT/lib/"libwallClockTimeStatistics_*.so > /dev/null
# Both arms run the case fvSolution, which must be the shipped fvSolution.fixedIter
cmp "$SHIP/system/fvSolution.fixedIter" "$CASE/system/fvSolution"
echo "== fvSolution: case = shipped fvSolution.fixedIter (cmp), both arms" >> compliance.txt
if [[ -d $CASE/postProcessing ]]; then mv "$CASE/postProcessing" "$ROOT/postProcessing.before"; fi

# --- solver environment ------------------------------------------------------------
set +eu
module purge
if [[ $COMPILER == gcc ]]; then
    module load gcc/14.2.0 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
    export CC=gcc CXX=g++ I_MPI_CC=gcc I_MPI_CXX=g++
    source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Gcc WM_CC=gcc WM_CXX=g++ \
        WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
else
    module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
    export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
    source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
        WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
fi
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
set -eu
[[ $WM_OPTIONS == "$OPT" ]]
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
ucxT=(-genv UCX_ZCOPY_THRESH inf -genv UCX_RNDV_THRESH inf)
{
    echo "simpleFoam=$(command -v simpleFoam)"
    echo "compiler=$($CXX --version | head -1) (COMPILER=$COMPILER, platform $OPT)"
    echo "mpi=$(mpirun --version 2>&1 | head -1)"
    echo "ucx=$(ucx_info -v | head -1) env: ${ucxT[*]}"
    echo "fixedIterSolvers=$FIS sha256=$(sha256sum "$FIS" | cut -d' ' -f1)"
    echo "fused=$FLIB sha256=$(sha256sum "$FLIB" | cut -d' ' -f1)"
    echo "stepReuse=$RLIB sha256=$(sha256sum "$RLIB" | cut -d' ' -f1)"
    echo "compactAddr=$CALIB sha256=$(sha256sum "$CALIB" | cut -d' ' -f1)"
    echo "simpleFoamFused=$SFF sha256=$(sha256sum "$SFF" | cut -d' ' -f1)"
} > build-metadata.txt
mpirun -np $NP -ppn 104 hostname | sort | uniq -c > placement.txt
if [[ $(wc -l < placement.txt) != $NNODES ]] || awk '$1 != 104 {bad=1} END {exit !bad}' placement.txt; then
    echo "PLACEMENT_WRONG:"; cat placement.txt; exit 1
fi

absent() {   # $1 text that must not be in log $2
    if grep -qF "$1" "$2"; then echo "UNEXPECTED in $2: $1"; return 1; fi
}
run() {   # $1 arm (basic|bonus), $2 repeat
    local a=$1 r=$2 dir=$ROOT/run_$1_$2 status=0
    mkdir "$dir"
    cd "$CASE"
    cp "$ROOT/controlDict.run" system/controlDict
    local pre=() app=simpleFoam
    if [[ $a == bonus ]]; then
        pre=(-genv LD_PRELOAD "$BONUS_PRELOAD" -genv STEP_REUSE_OFF "$BONUS_OFF" -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1 -genv SIMPLEFOAM_FUSED_REPORT 1)
        app=$SFF
    fi
    cp system/controlDict system/fvSolution "$dir/"
    echo "RUN_START arm=$a repeat=$r utc=$(date -u +%FT%TZ)"
    timeout --kill-after=30s 2400s mpirun -np $NP -ppn 104 "${glibcT[@]}" "${ucxT[@]}" "${pre[@]}" "$app" -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cp "$ROOT/controlDict.run" system/controlDict
    cd "$ROOT"
    # the run used the shipped fvSolution and the run controlDict (shipped controlDict.noWrite, endTime only)
    cmp "$SHIP/system/fvSolution.fixedIter" "$dir/fvSolution"
    cmp controlDict.run "$dir/controlDict"
    [[ $status == 0 ]]
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    if [[ $a == bonus ]]; then
        local log=$dir/simpleFoam.log f
        if [[ $FIXED_ITER == 1 ]]; then
            for f in k omega; do grep -qF "fixedIterSolvers: smoothSolver $f: fixed-iteration path" "$log"; done
        else
            for f in k omega; do grep -qF "fixedIterSolvers: smoothSolver $f: stock solve" "$log"; done
            absent 'fixed-iteration path' "$log"
        fi
        if [[ $FUSED_U == 1 ]]; then
            grep -qF 'fusedVectorSolve: U fused path (v3 folded copy) with every per-iteration residual check (v5)' "$log"
            [[ $(grep -c 'fusedVectorSolve: U fused path' "$log") == $(grep -c 'per-iteration residual check (v5)' "$log") ]]
        else
            absent 'fusedVectorSolve: U fused path' "$log"
        fi
        if [[ $GCACHE == 1 ]]; then
            for f in p k omega; do grep -qF "fixedIterSolvers: cache: grad($f) cached in code" "$log"; done
        else
            absent 'cached in code' "$log"
        fi
        # the folded Gauss-Seidel runs in compactAddr (which replaces stepReuse's worker) or in stepReuse
        grep -qE 'stepReuse: GaussSeidel: folded|compactAddr: GaussSeidel: 16-bit neighbour offsets, folded source copy' "$log"
        grep -q 'stepReuse: kOmegaSST: shared divU and k/omega convection, fused per-cell maths' "$log"
        grep -q 'stepReuse: divDevReff: fused explicit term' "$log"
        grep -q 'stepReuse: laplacian: fused Gauss linear assembly' "$log"
        grep -q 'stepReuse: convection: fused Gauss linearUpwindV' "$log"
        grep -q 'stepReuse: gaussGrad: interpolation fused' "$log"
        grep -q 'stepReuse: cellLimitedGrad: limiter passes fused' "$log"
        absent 'stepReuse: forceCoeffs: coefficient fields only with writeFields' "$log"
        absent 'stepReuse: wallShearStress: evaluated in write' "$log"
        grep -q 'simpleFoamFused: momentum solve without matrix copy' "$log"
        grep -q 'simpleFoamFused: pressure equation steps fused' "$log"
        grep -q 'compactAddr: GaussSeidel: 16-bit neighbour offsets' "$log"
        grep -q 'compactAddr: Amul: 16-bit neighbour offsets' "$log"
    fi
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS arm=$a repeat=$r Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
for ((r = 1; r <= REPEATS; r++)); do for a in $ARMS; do run $a $r; done; done
cmp "$CASE/system/controlDict" controlDict.run   # case restored

python3 - "$ROOT" "$NNODES" "$STEPS" "$REPEATS" "$ARMS" "$REF_CD" "GCACHE=$GCACHE FIXED_ITER=$FIXED_ITER FUSED_U=$FUSED_U STEP_REUSE_OFF=$BONUS_OFF" <<'PY'
import re, sys
from pathlib import Path
root, nn, steps, reps = Path(sys.argv[1]), sys.argv[2], sys.argv[3], int(sys.argv[4])
arms, ref_cd, switches = sys.argv[5].split(), sys.argv[6], sys.argv[7]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
num = re.compile(r'[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?')
def log(a, r): return (root / f'run_{a}_{r}' / 'simpleFoam.log').read_text(errors='replace')
def t(a, r): return float(pat.findall(log(a, r))[-1])
def rows(a, r, fo, name):
    return [l for p in sorted((root / f'run_{a}_{r}' / 'postProcessing' / fo).glob(f'*/{name}*.dat'))
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
def cd(a, r): return rows(a, r, 'forceCoeffsAll', 'coefficient')[-1].split()[1]
def residuals(a, r): return [l for l in log(a, r).splitlines() if 'Solving for' in l or 'continuity errors' in l]
def types(a, r):   # line types (numbers -> #) without the application line (simpleFoam / simpleFoamFused)
    return {' '.join(num.sub('#', l).split()) for l in log(a, r).splitlines() if l.strip() and not l.startswith('Exec ')}
both = 'basic' in arms and 'bonus' in arms
out = [f'# Final runs, {nn} node(s) x 104 ranks, occDrivAerStaticMesh 65M, {steps} steps, {reps} repeat(s), arms: {" ".join(arms)}', '',
       'basic = Basic tier (simpleFoam, shipped fvSolution and controlDict). bonus = simpleFoamFused + LD_PRELOAD compactAddr : fixedIterSolvers'
       + (' (fixed-iteration smoothSolver path, grad(p), grad(k), grad(omega) cached in code)'
          if 'FIXED_ITER=1' in switches else
          ' (grad(p), grad(k), grad(omega) cached in code; solver path off: k and omega on the stock smoothSolver)')
       + (' : fusedVectorSolve v5 (every per-iteration check computed)' if 'FUSED_U=1' in switches else '') + ' : stepReuse' +
       f', with the same shipped fvSolution and controlDict; bonus switches: {switches}. Arms interleaved on the same nodes.', '',
       '| Repeat | ' + ' | '.join(f'{a} s/step' for a in arms) + (' | bonus vs basic' if both else '') + ' | ' + ' | '.join(f'Cd {a} (step {steps})' for a in arms) + ' | Cd check |',
       '|---:|' + '---:|' * (len(arms) * 2 + (1 if both else 0)) + ':---|']
checks, ok = [], True
for r in range(1, reps + 1):
    ts = [t(a, r) for a in arms]; cds = [cd(a, r) for a in arms]
    ref = cd('basic', r) if both else ref_cd
    check = ('identical' if all(c == ref for c in cds) else 'DIFFERENT') + ('' if both else f' vs REF_CD {ref}' if ref else '')
    ok &= not ref or all(c == ref for c in cds)
    out.append(f'| {r} | ' + ' | '.join(f'{x:.4f}' for x in ts) + (f' | {100 * (t("bonus", r) / t("basic", r) - 1):+.2f}%' if both else '')
               + ' | ' + ' | '.join(cds) + f' | {check if ref else "n/a"} |')
    if both:
        same = {'coefficient.dat': rows('basic', r, 'forceCoeffsAll', 'coefficient') == rows('bonus', r, 'forceCoeffsAll', 'coefficient'),
                'solverInfo.dat': rows('basic', r, 'solverInfo', 'solverInfo') == rows('bonus', r, 'solverInfo', 'solverInfo'),
                'residual lines': residuals('basic', r) == residuals('bonus', r)}
        missing = sorted(types('basic', r) - types('bonus', r))
        ok &= all(same.values()) and not missing
        checks.append(f'repeat {r}: ' + ', '.join(f'{k} {"identical" if v else "DIFFERENT"}' for k, v in same.items())
                      + f' ({len(rows("basic", r, "forceCoeffsAll", "coefficient"))} coefficient rows)')
        checks.append(f'repeat {r}: F5 ' + ('every basic log line type present in the bonus log' if not missing
                                            else 'basic line types MISSING in the bonus log: ' + ' | '.join(missing)))
best = {a: min(t(a, r) for r in range(1, reps + 1)) for a in arms}
out += [''] + [f'Best {a}: **{v:.4f} s/step**' for a, v in best.items()]
if checks: out += ['', 'Bonus vs basic, same job:', ''] + [f'- {c}' for c in checks]
(root / 'REPORT.md').write_text('\n'.join(out) + '\n'); print('\n'.join(out))
with open(root / 'compliance.txt', 'a') as f:
    f.write('== bonus vs basic (F3 self-proof, F5)\n' + ''.join(f'   {c}\n' for c in checks))
sys.exit(0 if ok else 1)
PY
echo "FINAL_COMPLETE root=$ROOT"
