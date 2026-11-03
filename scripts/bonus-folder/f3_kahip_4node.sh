#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N f3_kahip_4n
#PBS -j oe

# F3 self-proof for the bonus folder's KaHIP configuration at 4 nodes x 104 ranks, 150 steps (Addendum F3: the
# baseline is "the same configuration as your Bonus run (same decomposition method and coefficients, same rank count,
# same node count), with the change removed"). Both arms on the same linked copy of the team KaHIP partition
# (an0642 kahip_4_node via kahip/regen_bench_181335722/case_kahip, decomposed with the bonus folder's multiLevel
# strong/eco dictionary), the bonus folder's dicts/4node, and the solver environment of its run/4node/run.sh
# (bonus-folder build, -xCORE-AVX512; malloc and UCX settings; -fileHandler hostCollated; mpirun -np 416):
#   base  : stock simpleFoam of that build (the change removed)
#   bonus : simpleFoamFused + LD_PRELOAD compactAddr:fixedIterSolvers:fusedVectorSolve:stepReuse, STEP_REUSE_OFF=fixed
#           (run.sh since 2026-11-03: every per-iteration convergence check computed: fusedVectorSolve v5 for U, k and
#           omega on the stock smoothSolver; job 182978547 and earlier ran fusedVectorSolve v3 and the fixed-iteration
#           path, which skip those checks, with STEP_REUSE_OFF=wss,fco)
# Pass: coefficient.dat bit-identical (else every Cd within +-0.5%). Also solverInfo.dat and residual lines, F5 (every
# base log line type in the bonus log), and the Cd sanity range (150 steps: 0.284-0.292, "about 0.28-0.29").
# 1 and 2 nodes are not run: no KaHIP partitions exist for them and the user chose not to decompose.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
BON=/scratch/ya41/apac-2026-bonus/openfoam
RES=/scratch/ya41/kn3896/openfoam/contest-results
REF=$RES/compliance/reference/run/case/occDrivAerStaticMesh/system
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
SCOTCH=$RES/decomp-screen/case416_w05
KAHIP=$RES/kahip/regen_bench_181335722/case_kahip
ROOT=$RES/bonus-folder/f3_kahip_4node_${PBS_JOBID%%.*}
STEPS=150

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "F3_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$BON/run/4node/intel-scatter-cpus.awk" "$BON/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
cp "$BON/run/4node/run.sh" run.sh.bonus-folder
( cd "$BON" && git log -1 --format='bonus folder repo HEAD %h %s' && git status --short ) > bonus-folder-state.txt 2>&1 || true

# --- case copies -------------------------------------------------------------------------
mkcase() {   # $1 dir, $2 prepared case (processors*), $3 caseDefinition, $4 decomposeParDict, $5 controlDict
    local d=$1
    mkdir -p "$d/constant" "$d/system"
    cp -a "$BON/run/case/occDrivAerStaticMesh/system/." "$d/system/"
    cp "$3" "$d/system/include/caseDefinition"
    cp "$4" "$d/system/decomposeParDict"
    cp "$5" "$d/system/controlDict"
    cp -a "$BON"/run/case/occDrivAerStaticMesh/constant/{transportProperties,turbulenceProperties} "$d/constant/"
    ln -s "$MESH" "$d/constant/polyMesh"
    for p in "$2"/processors416_*; do ln -s "$p" "$d/${p##*/}"; done
}
mkcase "$ROOT/case_kahip" "$KAHIP" "$BON/dicts/4node/caseDefinition" "$BON/dicts/4node/decomposeParDict" \
    "$BON/dicts/4node/controlDict"
# Compliance of the case files against the shipped case (Addendum F1; Bonus may change the method)
{
    for a in kahip; do
        s=$ROOT/case_$a/system
        cmp -s "$REF/fvSolution.fixedIter" "$s/fvSolution" || echo "$a: fvSolution DIFFERS from shipped fvSolution.fixedIter"
        cmp -s "$REF/fvSchemes" "$s/fvSchemes" || echo "$a: fvSchemes DIFFERS"
        echo "== $a controlDict vs shipped controlDict.noWrite:"; diff "$REF/controlDict.noWrite" "$s/controlDict" || true
        echo "== $a caseDefinition vs shipped:"; diff "$REF/include/caseDefinition" "$s/include/caseDefinition" || true
    done
} > compliance.txt
cat compliance.txt
[[ $(diff "$REF/controlDict.noWrite" "$ROOT/case_kahip/system/controlDict" | grep -c '^[<>]') == 2 ]]
grep -q "^endTime *$STEPS;" "$ROOT/case_kahip/system/controlDict"

# --- solver environment: run.sh lines after the preparation chain ------------------------------
set +eu
module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2
module load metis/5.1.0 fftw3/3.3.10 boost/1.84.0
source "$BON/src/OpenFOAM-v2512/etc/bashrc"
export WM_COMPILER_TYPE=system
export WM_COMPILER=Icx
export WM_MPLIB=INTELMPI
export WM_COMPILE_OPTION=Opt
export WM_PRECISION_OPTION=DP
export WM_LABEL_SIZE=32
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export FOAM_SIGFPE=true
export I_MPI_LIBRARY_KIND=release
export I_MPI_FABRICS=shm:ofi
export I_MPI_OFI_PROVIDER=mlx
export I_MPI_DEBUG=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_EXCLUDE_LIST I_MPI_PIN_DOMAIN_SIZE
export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core
export I_MPI_HYDRA_TOPOLIB=hwloc
I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
export I_MPI_PIN_PROCESSOR_LIST
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-4n-104ppn-scatter.dat
export UCX_ZCOPY_THRESH=inf UCX_RNDV_THRESH=inf
export MALLOC_MMAP_THRESHOLD_=33554432 MALLOC_TRIM_THRESHOLD_=4294967296 MALLOC_TOP_PAD_=268435456   # run.sh line 84 since 2026-11-02
REPO=$(cd "$WM_PROJECT_DIR/.." && pwd)/apac-hpc-openfoam   # as build.pbs/run.sh since 2026-11-03 (logical path)
LIB=$REPO/platforms/$WM_OPTIONS/lib
PRELOAD=$LIB/libcompactAddr.so:$LIB/libfixedIterSolvers.so:$LIB/libfusedVectorSolve.so:$LIB/libstepReuse.so   # as run.sh since 2026-11-03
SOLVER=$REPO/platforms/$WM_OPTIONS/bin/simpleFoamFused
parEx="mpirun -np 416"
fileHandler="-fileHandler hostCollated"
set -eu
{
    echo "WM_PROJECT_DIR=$WM_PROJECT_DIR -> $(readlink -f "$WM_PROJECT_DIR")"
    echo "WM_OPTIONS=$WM_OPTIONS"
    echo "simpleFoam=$(command -v simpleFoam) sha256=$(sha256sum "$(command -v simpleFoam)" | cut -c1-16)"
    echo "simpleFoamFused=$SOLVER sha256=$(sha256sum "$SOLVER" | cut -c1-16)"
    for l in ${PRELOAD//:/ }; do echo "preload $l sha256=$(sha256sum "$l" | cut -c1-16)"; done
    echo "clone HEAD=$(git -C "$REPO" rev-parse --short HEAD)"
    echo "cxxflags=$(wmake -show-cxxflags 2>/dev/null)"
    echo "mpi=$(mpirun --version 2>&1 | head -1)"
} > build-metadata.txt
cat build-metadata.txt
[[ $WM_OPTIONS == linux64IcxDPInt32Opt && -x $SOLVER && -x $(command -v simpleFoam) ]]

run() {   # $1 arm
    local a=$1 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    cd "$ROOT/case_kahip"
    echo "RUN_START arm=$a utc=$(date -u +%FT%TZ)"
    if [[ $a == bonus ]]; then
        STEP_REUSE_OFF=fixed STEP_REUSE_REPORT=1 FUSED_VECTOR_SOLVE_REPORT=1 timeout --kill-after=30s 1500s \
            $parEx -genv LD_PRELOAD "$PRELOAD" -genv STEP_REUSE_OFF fixed "$SOLVER" -parallel $fileHandler > "$dir/simpleFoam.log" 2>&1 || status=$?
    else
        timeout --kill-after=30s 1500s $parEx simpleFoam -parallel $fileHandler > "$dir/simpleFoam.log" 2>&1 || status=$?
    fi
    [[ -d postProcessing ]] && mv postProcessing "$dir/"
    cd "$ROOT"
    [[ -z $(ls -d "$ROOT/case_kahip"/processors416_*/[1-9]* 2>/dev/null) ]]   # nothing written into the prepared case
    [[ $status == 0 ]] || { echo "RUN_FAILED arm=$a status=$status"; tail -20 "$dir/simpleFoam.log"; return 1; }
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
run base
run bonus

# --- checks --------------------------------------------------------------------------------------
L=$ROOT/run_bonus/simpleFoam.log
echo "== bonus report lines (master):"
grep -hE '^(fixedIterSolvers|stepReuse|compactAddr|fusedVectorSolve|simpleFoamFused):' "$L" | sort | uniq -c
miss=0
for p in 'fixedIterSolvers: smoothSolver k: stock solve' 'fixedIterSolvers: smoothSolver omega: stock solve' \
         'fusedVectorSolve: U fused path (v3 folded copy) with every per-iteration residual check (v5)' \
         'fixedIterSolvers: cache: grad(p) cached in code' 'fixedIterSolvers: cache: grad(k) cached in code' \
         'fixedIterSolvers: cache: grad(omega) cached in code' \
         'compactAddr: GaussSeidel: 16-bit neighbour offsets' 'compactAddr: Amul: 16-bit neighbour offsets' \
         'stepReuse: kOmegaSST' 'stepReuse: divDevReff: fused' 'stepReuse: laplacian: fused' 'stepReuse: convection: fused' \
         'stepReuse: gaussGrad: interpolation fused' 'stepReuse: cellLimitedGrad'; do
    grep -qF "$p" "$L" || { echo "REPORT MISSING: $p"; miss=1; }
done
for p in 'stepReuse: wallShearStress: evaluated in write' 'stepReuse: forceCoeffs: coefficient fields only' \
         'fixed-iteration path' 'smoothSolver Ux: stock solve'; do
    grep -qF "$p" "$L" && { echo "REPORT UNEXPECTED: $p"; miss=1; }
done
python3 - "$ROOT" "$miss" <<'PY'
import re, sys, glob
from pathlib import Path
root, miss = Path(sys.argv[1]), sys.argv[2] == '1'
def rows(a, obj, f):
    p = glob.glob(str(root / f'run_{a}/postProcessing/{obj}/*/{f}'))[0]
    return [l for l in open(p) if l.strip() and not l.startswith('#')]
def t(a):
    return float(re.findall(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)', (root / f'run_{a}/simpleFoam.log').read_text())[-1])
def res(a):
    return [l for l in open(root / f'run_{a}/simpleFoam.log', errors='replace') if re.search(r'Solving for|continuity errors', l)]
num = re.compile(r'[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?')
def types(a):
    s = set()
    for l in open(root / f'run_{a}/simpleFoam.log', errors='replace'):
        l = ' '.join(num.sub('#', l).split())
        if l and not l.startswith('Exec ') and not l.startswith('Case '):
            s.add(l)
    return s
cb, cn = rows('base', 'forceCoeffsAll', 'coefficient*.dat'), rows('bonus', 'forceCoeffsAll', 'coefficient*.dat')
sb, sn = rows('base', 'solverInfo', 'solverInfo.dat'), rows('bonus', 'solverInfo', 'solverInfo.dat')
cd_b, cd_n = float(cb[-1].split()[1]), float(cn[-1].split()[1])
same = cb == cn
worst = max(abs(float(n.split()[1]) / float(b.split()[1]) - 1) * 100 for b, n in zip(cb, cn) if float(b.split()[1]) != 0)
missing = sorted(types('base') - types('bonus'))
tb, tn = t('base'), t('bonus')
ok = same or abs(cd_n / cd_b - 1) * 100 <= 0.5
out = [
 '# F3 self-proof, bonus folder KaHIP configuration, 4 nodes x 104, occDrivAerStaticMesh 65M, 150 steps, same job', '',
 '| Arm | s/step | Cd step 150 |', '|---|---:|---|',
 f'| base: stock simpleFoam, team KaHIP partition | {tb:.4f} | {cd_b:.17g} |',
 f'| bonus: bonus folder run.sh configuration, same partition | {tn:.4f} | {cd_n:.17g} |', '',
 f'F3 (coefficient.dat bit-identical, else Cd within +-0.5%): ' + ('PASS' if ok else 'FAIL') +
   f' -- coefficient.dat identical: {same} ({len(cb)} rows; largest Cd deviation {worst:.2e}%)',
 f'solverInfo.dat identical: {sb == sn}; residual lines identical: {res("base") == res("bonus")}',
 'F5 (every base log line type in the bonus log): ' + ('PASS' if not missing else f'FAIL ({len(missing)} missing)'),
 'report lines (stepReuse, fixedIterSolvers, compactAddr): ' + ('PASS' if not miss else 'see job.log'),
 f'Sanity check, Cd at 150 steps {cd_n:.4f} vs references 0.284-0.292 ("about 0.28-0.29")',
 f'Bonus vs base (same partition): {100*(tn/tb-1):+.2f}%; vs the Basic result of this build on scotch (1.9360, job 182944449): {100*(tn/1.9360-1):+.2f}%',
]
if missing:
    out += ['missing line types:'] + ['    ' + m for m in missing[:40]]
(root / 'REPORT.md').write_text('\n'.join(out) + '\n'); print('\n'.join(out))
PY
python3 /scratch/ya41/kn3896/openfoam/tools/perfModel/ab_steps.py "$ROOT" 2>/dev/null | tail -6 || true
# pack the case copies (inode quota): config and dynamicCode only, links not followed
for d in case_kahip; do
    n=$(find "$d" | wc -l); tar czf "$d.tar.gz" "$d" && [[ $(tar tzf "$d.tar.gz" | wc -l) == "$n" ]] && rm -rf "$d"
done
echo "F3_COMPLETE root=$ROOT"
