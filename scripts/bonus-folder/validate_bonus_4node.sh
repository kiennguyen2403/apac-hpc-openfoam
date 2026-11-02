#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N validate_bonus_4n
#PBS -j oe

# Validation of the bonus folder (/scratch/ya41/apac-2026-bonus/openfoam) at 4 nodes x 104 ranks, 150 steps.
# Solver environment exactly as its run/4node/run.sh after the preparation steps: its modules, its etc/bashrc
# (../../src/OpenFOAM-v2512 -> Basic folder's tree, built by build/build.pbs with prefs.sh -xCORE-AVX512), its
# Intel MPI pinning and tuned Allreduce files, UCX setting, STEP_REUSE_OFF=wss,fco, -fileHandler hostCollated,
# "mpirun -np 416". Two arms, same job, fixed order:
#   basic : stock simpleFoam of that build, scotch (4 104)/(1 0.5) partition (prepared case decomp-screen/case416_w05)
#           = the F3 baseline (Basic tier allows only scotch/hierarchical)
#   bonus : simpleFoamFused + LD_PRELOAD compactAddr:fixedIterSolvers:fusedVectorSolve:stepReuse from the clone
#           ../src/apac-hpc-openfoam (bf45d4c), its dicts/4node (multiLevel KaHIP strong/eco), on the team's KaHIP
#           partition (an0642 kahip_4_node, decomposed with the same dictionary, hostCollated, renumbered as in
#           Allrun; fields regenerated in kahip/regen_bench_181335722/case_kahip).
# Deviation from run.sh: the preparation chain is not repeated (the KaHIP decomposition alone takes 1.3-1.7 h serial,
# job 182167809, more than run.sh's 40-minute walltime, and the bonus case has no mesh; its tarball is gone).
# Both arms run on linked copies (processors* symlinked, controlDict noWrite: nothing is written into them).
# Checks: F3 (Cd at step 150 within +-0.5% of basic; whole coefficient.dat compared), F5 (every basic log line type
# in the bonus log), report lines of every library, timing against the compliant finals (182882563).

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
BON=/scratch/ya41/apac-2026-bonus/openfoam
RES=/scratch/ya41/kn3896/openfoam/contest-results
REF=$RES/compliance/reference/run/case/occDrivAerStaticMesh/system
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
SCOTCH=$RES/decomp-screen/case416_w05
KAHIP=$RES/kahip/regen_bench_181335722/case_kahip
ROOT=$RES/bonus-folder/validate_4node_${PBS_JOBID%%.*}
STEPS=150

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "VALIDATE_FAILED line=$LINENO status=$?"' ERR
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
mkcase "$ROOT/case_basic" "$SCOTCH" "$SCOTCH/system/include/caseDefinition" "$SCOTCH/system/decomposeParDict" \
    "$BON/dicts/4node/controlDict"
mkcase "$ROOT/case_bonus" "$KAHIP" "$BON/dicts/4node/caseDefinition" "$BON/dicts/4node/decomposeParDict" \
    "$BON/dicts/4node/controlDict"
# Compliance of the case files against the shipped case (Addendum F1; Bonus may change the method)
{
    for a in basic bonus; do
        s=$ROOT/case_$a/system
        cmp -s "$REF/fvSolution.fixedIter" "$s/fvSolution" || echo "$a: fvSolution DIFFERS from shipped fvSolution.fixedIter"
        cmp -s "$REF/fvSchemes" "$s/fvSchemes" || echo "$a: fvSchemes DIFFERS"
        echo "== $a controlDict vs shipped controlDict.noWrite:"; diff "$REF/controlDict.noWrite" "$s/controlDict" || true
        echo "== $a caseDefinition vs shipped:"; diff "$REF/include/caseDefinition" "$s/include/caseDefinition" || true
    done
} > compliance.txt
cat compliance.txt
[[ $(diff "$REF/controlDict.noWrite" "$ROOT/case_bonus/system/controlDict" | grep -c '^[<>]') == 2 ]]
grep -q "^endTime *$STEPS;" "$ROOT/case_bonus/system/controlDict"

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
REPO=$WM_PROJECT_DIR/../apac-hpc-openfoam
LIB=$REPO/platforms/$WM_OPTIONS/lib
PRELOAD=$LIB/libcompactAddr.so:$LIB/libfixedIterSolvers.so:$LIB/libfusedVectorSolve.so:$LIB/libstepReuse.so
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
    cd "$ROOT/case_$a"
    echo "RUN_START arm=$a utc=$(date -u +%FT%TZ)"
    if [[ $a == bonus ]]; then
        STEP_REUSE_OFF=wss,fco STEP_REUSE_REPORT=1 timeout --kill-after=30s 1500s \
            $parEx -genv LD_PRELOAD "$PRELOAD" "$SOLVER" -parallel $fileHandler > "$dir/simpleFoam.log" 2>&1 || status=$?
    else
        timeout --kill-after=30s 1500s $parEx simpleFoam -parallel $fileHandler > "$dir/simpleFoam.log" 2>&1 || status=$?
    fi
    [[ -d postProcessing ]] && mv postProcessing "$dir/"
    cd "$ROOT"
    [[ -z $(ls -d "$ROOT/case_$a"/processors416_*/[1-9]* 2>/dev/null) ]]   # nothing written into the prepared cases
    [[ $status == 0 ]] || { echo "RUN_FAILED arm=$a status=$status"; tail -20 "$dir/simpleFoam.log"; return 1; }
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $STEPS ]]
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS arm=$a Cd_step$STEPS=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
run basic
run bonus

# --- checks --------------------------------------------------------------------------------------
L=$ROOT/run_bonus/simpleFoam.log
echo "== bonus report lines (master):"
grep -hE '^(fixedIterSolvers|stepReuse|compactAddr|fusedVectorSolve|simpleFoamFused):' "$L" | sort | uniq -c
miss=0
for p in 'fixedIterSolvers: smoothSolver k: fixed-iteration path' 'fixedIterSolvers: smoothSolver omega: fixed-iteration path' \
         'fixedIterSolvers: cache: grad(p) cached in code' 'fixedIterSolvers: cache: grad(k) cached in code' \
         'fixedIterSolvers: cache: grad(omega) cached in code' 'fusedVectorSolve: U fused path' \
         'compactAddr: GaussSeidel: 16-bit neighbour offsets' 'compactAddr: Amul: 16-bit neighbour offsets' \
         'stepReuse: kOmegaSST' 'stepReuse: divDevReff: fused' 'stepReuse: laplacian: fused' 'stepReuse: convection: fused' \
         'stepReuse: gaussGrad: interpolation fused' 'stepReuse: cellLimitedGrad' 'simpleFoamFused: momentum solve' \
         'simpleFoamFused: pressure equation steps fused'; do
    grep -qF "$p" "$L" || { echo "REPORT MISSING: $p"; miss=1; }
done
for p in 'stepReuse: wallShearStress: evaluated in write' 'stepReuse: forceCoeffs: coefficient fields only' 'stock solve'; do
    grep -qF "$p" "$L" && { echo "REPORT UNEXPECTED: $p"; miss=1; }
done
python3 - "$ROOT" "$miss" <<'PY'
import re, sys, glob
from pathlib import Path
root, miss = Path(sys.argv[1]), sys.argv[2] == '1'
def coeffs(a):
    f = glob.glob(str(root / f'run_{a}/postProcessing/forceCoeffsAll/*/coefficient*.dat'))[0]
    return [l.split() for l in open(f) if l.strip() and not l.startswith('#')]
def t(a):
    return float(re.findall(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)', (root / f'run_{a}/simpleFoam.log').read_text())[-1])
num = re.compile(r'[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?')
def types(a):
    s = set()
    for l in open(root / f'run_{a}/simpleFoam.log', errors='replace'):
        l = ' '.join(num.sub('#', l.replace(str(root / f'case_{a}'), 'CASE')).split())
        if l and not l.startswith('Exec ') and not l.startswith('Case '):
            s.add(l)
    return s
cb, cn = coeffs('basic'), coeffs('bonus')
cd_b, cd_n = float(cb[-1][1]), float(cn[-1][1])
rel = 100 * (cd_n / cd_b - 1)
worst = max(abs(float(n[1]) / float(b[1]) - 1) * 100 for b, n in zip(cb, cn) if float(b[1]) != 0)
missing = sorted(types('basic') - types('bonus'))
tb, tn = t('basic'), t('bonus')
out = [
 '# Bonus folder validation, 4 nodes x 104, occDrivAerStaticMesh 65M, 150 steps, same job', '',
 '| Arm | s/step | vs basic | compliant finals 182882563 | Cd step 150 |', '|---|---:|---:|---:|---|',
 f'| basic (stock simpleFoam, this build, scotch) | {tb:.4f} | - | 1.9536 | {cd_b:.17g} |',
 f'| bonus (bonus folder run.sh config, KaHIP) | {tn:.4f} | {100*(tn/tb-1):+.2f}% | 1.2674 | {cd_n:.17g} |', '',
 f'F3: Cd at step 150 bonus vs basic {rel:+.4f}% (limit +-0.5%): ' + ('PASS' if abs(rel) <= 0.5 else 'FAIL'),
 f'    largest Cd deviation over all 150 rows: {worst:.4f}%; coefficient.dat identical: {cb == cn}',
 'F5: every basic log line type present in the bonus log: ' + ('PASS' if not missing else f'FAIL ({len(missing)} missing)'),
 'report lines of every library: ' + ('PASS' if not miss else 'FAIL (see job.log)'),
]
if missing:
    out += ['missing line types:'] + ['    ' + m for m in missing[:40]]
(root / 'REPORT.md').write_text('\n'.join(out) + '\n'); print('\n'.join(out))
PY
python3 /scratch/ya41/kn3896/openfoam/tools/perfModel/ab_steps.py "$ROOT" 2>/dev/null | tail -6 || true
# pack the case copies (inode quota): config and dynamicCode only, links not followed
for d in case_basic case_bonus; do
    n=$(find "$d" | wc -l); tar czf "$d.tar.gz" "$d" && [[ $(tar tzf "$d.tar.gz" | wc -l) == "$n" ]] && rm -rf "$d"
done
echo "VALIDATE_COMPLETE root=$ROOT"
