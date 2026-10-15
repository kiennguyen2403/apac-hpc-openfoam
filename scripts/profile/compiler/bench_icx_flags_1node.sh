#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=01:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N icx_flags_1n
#PBS -j oe

# Phase 1: ICX flag tuning on occDrivAerStaticMesh 65M, 1 node, 104 ranks, 40 steps.
# Variants (custom v2512 + Intel MPI 2021.17.2, all -xSAPPHIRERAPIDS -fp-model=precise):
#   base = OptSPR (-O3)   zmm = OptSPRzmm (-O3 -qopt-zmm-usage=high)
#   O2 = OptSPRO2 (-O2)   ipo = OptSPRipo (-O3 -ipo -ffat-lto-objects)
# Reuses the prepared case (scotch 104, hostCollated, renumbered, potentialFoam,
# applyBoundaryLayer) -- no decomposition or prep here. Team runtime settings:
# scatter pinning + intelmpi-allreduce-1n-104ppn-scatter.dat.
# Per variant: 3-step warmup (compiles the coded timer), then 2 timed repeats, mirrored order.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
NP=104
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
TEAM=/scratch/ya41/apac-2026/openfoam
CASE=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/icx-flags/drivaer_1node_${PBS_JOBID%%.*}
VARIANTS=(base zmm O2 ipo)
fh="-fileHandler hostCollated"

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "PHASE1_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
cp "$TEAM/run/1node/intel-scatter-cpus.awk" "$TEAM/run/1node/intelmpi-allreduce-1n-104ppn-scatter.dat" .

# Reused prepared case: prep logs must be complete and endTime must be the 1-node 40.
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do
    grep -q '^End' "$CASE/../logFiles/$l.log"
done
cmp "$TEAM/dicts/1node/controlDict" "$CASE/system/controlDict"
[[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' "$CASE/system/controlDict") == 40 ]]
echo "REUSING prepared case $CASE"
[[ -d $CASE/postProcessing ]] && mv "$CASE/postProcessing" "$ROOT/postProcessing.before-${PBS_JOBID%%.*}"

load_variant() {   # $1 = base|zmm|O2|ipo ; runs in a subshell
    local tune=$1 opt=OptSPR
    [[ $tune != base ]] && opt=OptSPR$tune
    set +eu
    module purge
    module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
    export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
    source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
        WM_MPLIB=INTELMPI WM_COMPILE_OPTION=$opt
    export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
    set -eu
    [[ $WM_OPTIONS == linux64IcxDPInt32$opt ]] || { echo "ERROR: WM_OPTIONS=$WM_OPTIONS"; return 1; }
    export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FOAM_SIGFPE=true
    export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx I_MPI_DEBUG=0
    unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_EXCLUDE_LIST I_MPI_PIN_DOMAIN_SIZE
    export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
    I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
    export I_MPI_PIN_PROCESSOR_LIST
    unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
    export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-1n-104ppn-scatter.dat
}

# Skip variants whose build did not finish (e.g. a failed ipo build).
AVAILABLE=()
for v in "${VARIANTS[@]}"; do
    opt=OptSPR; [[ $v != base ]] && opt=OptSPR$v
    log=$(ls -t /scratch/ya41/kn3896/openfoam/logs/build/v2512-icx-spr$([[ $v != base ]] && echo "-$v")-intelmpi/build-*.log 2>/dev/null | head -1 || true)
    if [[ -x $FOAM/platforms/linux64IcxDPInt32$opt/bin/simpleFoam ]] && grep -q '^BUILD_SUCCESS' "$log" 2>/dev/null; then
        AVAILABLE+=("$v")
    else
        echo "SKIP variant=$v (build missing or not successful: ${log:-no log})"
    fi
done
echo "VARIANTS ${AVAILABLE[*]}"

# Every check exits explicitly: run() is called from "if", where set -e is suspended,
# so one failing variant is recorded and skipped instead of stopping the job (job 180741885).
run() {   # $1 variant, $2 label (warmup|t1|t2)
    local v=$1 label=$2 dir=$ROOT/run_${1}_$2 steps=40
    mkdir "$dir" || return 1
    (
        load_variant "$v" || exit 2
        cd "$CASE" || exit 2
        if [[ $label == warmup ]]; then
            steps=3
            {
                echo "variant=$v wm_options=$WM_OPTIONS"
                echo "simpleFoam=$(command -v simpleFoam)"
                echo "flags=$(wmake -show-cxxflags)"
                readelf -p .comment "$(command -v simpleFoam)" | grep -v '^$'
            } > "$dir/build-metadata.txt"
            ldd "$(command -v simpleFoam)" > "$dir/simpleFoam-ldd.txt" || exit 3
            if grep -q 'not found' "$dir/simpleFoam-ldd.txt"; then exit 3; fi
            # OpenFOAM only compiles the coded timer when its source SHA1 changes, not per platform:
            # with the source already present (from another variant) it just fails to load this
            # variant's lib and runs without the timer (job 180744480). Build it here explicitly.
            timer_dir=$CASE/dynamicCode/wallClockTimeStatistics
            timer_lib=$(ls "$timer_dir"/../platforms/$WM_OPTIONS/lib/libwallClockTimeStatistics_*.so 2>/dev/null || true)
            if [[ -d $timer_dir && -z $timer_lib ]]; then
                (cd "$timer_dir" && wmake libso) > "$dir/timer-build.log" 2>&1 || exit 3
                ls "$timer_dir"/../platforms/$WM_OPTIONS/lib/libwallClockTimeStatistics_*.so > /dev/null || exit 3
            fi
        fi
        # Fresh team controlDict every run; only the warmup shortens endTime (foamDictionary would reformat the file).
        cp "$TEAM/dicts/1node/controlDict" system/controlDict || exit 4
        [[ $steps == 40 ]] || sed -i "s/^endTime .*/endTime           $steps;/" system/controlDict
        [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == $steps ]] || exit 4
        echo "RUN_START variant=$v run=$label utc=$(date -u +%FT%TZ)"
        status=0
        timeout --kill-after=30s 1200s mpirun -np $NP simpleFoam -parallel $fh > "$dir/simpleFoam.log" 2>&1 || status=$?
        cp "$TEAM/dicts/1node/controlDict" system/controlDict
        if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
        if [[ $status != 0 ]]; then
            echo "RUN_FAILED variant=$v run=$label mpirun_status=$status: $(grep -m1 -E 'FOAM FATAL|sigFpe|Segmentation' "$dir/simpleFoam.log" | cut -c1-150)"
            exit 5
        fi
        grep -qx End "$dir/simpleFoam.log" || exit 6
        if ! grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"; then
            echo "RUN_FAILED variant=$v run=$label: contest timer missing"; exit 7
        fi
        [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $steps ]] || exit 6
        cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
        echo "RUN_SUCCESS variant=$v run=$label Cd_step$steps=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
    )
}

OK=()
for v in "${AVAILABLE[@]}"; do
    if run "$v" warmup; then OK+=("$v"); else echo "VARIANT_DROPPED variant=$v (warmup failed)"; fi
done
(( ${#OK[@]} > 0 ))
for v in "${OK[@]}"; do run "$v" t1 || echo "RUN_FAILED_RECORDED variant=$v run=t1"; done
for ((i=${#OK[@]}-1; i>=0; i--)); do run "${OK[$i]}" t2 || echo "RUN_FAILED_RECORDED variant=${OK[$i]} run=t2"; done
cmp "$TEAM/dicts/1node/controlDict" "$CASE/system/controlDict"   # case left as found

# --- summary --------------------------------------------------------------------
python3 - "$ROOT" "${OK[@]}" <<'EOF'
import re, sys
from pathlib import Path
root, variants = Path(sys.argv[1]), sys.argv[2:]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
def t(v, r): return float(pat.findall((root / f'run_{v}_{r}' / 'simpleFoam.log').read_text())[-1])
def cd(v, r):
    rows = [l.split() for p in (root / f'run_{v}_{r}' / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat')
            for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    return float(rows[-1][1])
res = {}
for v in variants:
    try:
        res[v] = ((t(v, 't1') + t(v, 't2')) / 2, t(v, 't1'), t(v, 't2'), cd(v, 't1'))
    except (FileNotFoundError, IndexError) as e:
        print('incomplete variant', v, e)
base = res['base'][0] if 'base' in res else None
lines = ['# ICX flag tuning, occDrivAerStaticMesh 65M, 1 node, 104 ranks, 40 steps', '',
         'All: -xSAPPHIRERAPIDS -fp-model=precise, Intel MPI 2021.17.2, scatter pinning, tuned Allreduce.', '',
         '| Variant | Repeat 1 (s) | Repeat 2 (s) | Mean (s/step) | vs base | Cd step 40 (0.43-0.44) |',
         '|---|---:|---:|---:|---:|---:|']
for v in sorted(res, key=lambda k: res[k][0]):
    m, a, b, c = res[v]
    rel = '%+.2f%%' % (100 * (m / base - 1)) if base else 'n/a'
    ok = 'OK' if 0.43 <= c <= 0.44 else 'OUT OF RANGE'
    lines.append(f'| {v} | {a:.4f} | {b:.4f} | {m:.4f} | {rel} | {c:.4f} {ok} |')
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "PHASE1_COMPLETE root=$ROOT"
