#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N alloc_thp_1n
#PBS -j oe

# Memory allocator / transparent huge page test, occDrivAerStaticMesh 65M, 1 node, 104 ranks,
# 40 steps, ICX OptSPR + Intel MPI with the team runtime settings. Gadi THP mode is "madvise",
# so glibc malloc gets 4 KiB pages only. Variants (one run each, no mirror):
#   glibc   : default (reference)
#   glibcT  : glibc without mmap/munmap churn of large temporaries (MALLOC_* env)
#   jemalloc: LD_PRELOAD jemalloc 5.3.0, thp:always (madvise huge pages)
#   mimalloc: LD_PRELOAD mimalloc 3.1.5, allow_large_os_pages + eager arena commit
#   tbbmalloc: LD_PRELOAD Intel oneTBB 2022.3.1 tbbmalloc_proxy (module intel-tbb/2022.3.1)
# The first variant listed is the reference in REPORT.md.
# Reuses the prepared 104-rank case (no decomposition/prep). Before the runs, a probe checks
# how many huge pages this compute node grants; during each run the ranks' AnonHugePages is
# sampled from /proc so the effect can be tied to actual huge-page use.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
NP=104
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
TEAM=/scratch/ya41/apac-2026/openfoam
CASE=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case
T=/scratch/ya41/kn3896/openfoam/tools
JE=$T/jemalloc-5.3.0/lib/libjemalloc.so.2
MI=$T/mimalloc-3.1.5/lib64/libmimalloc.so.3
TBB=/apps/intel-tools/intel-tbb/2022.3.1/lib/libtbbmalloc_proxy.so.2
PROBE=/scratch/ya41/kn3896/openfoam/tools/thp-probe
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/alloc-thp/drivaer_1node_${PBS_JOBID%%.*}
fh="-fileHandler hostCollated"

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "ALLOC_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
cp "$TEAM/run/1node/intel-scatter-cpus.awk" "$TEAM/run/1node/intelmpi-allreduce-1n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/../logFiles/$l.log"; done
cp "$TEAM/dicts/1node/controlDict" "$CASE/system/controlDict"
[[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' "$CASE/system/controlDict") == 40 ]]
ls "$CASE/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
if [[ -d $CASE/postProcessing ]]; then mv "$CASE/postProcessing" "$ROOT/postProcessing.before"; fi
echo "REUSING prepared case $CASE"

# --- probe: huge pages this node grants (1 GiB, direct madvise and per allocator) --------
{
    cat /sys/kernel/mm/transparent_hugepage/enabled /sys/kernel/mm/transparent_hugepage/defrag
    grep -E 'MemFree|AnonHugePages' /proc/meminfo
    echo "direct madvise:   $("$PROBE/direct" | grep AnonHuge)"
    echo "glibc:            $("$PROBE/thp" | grep AnonHuge)"
    echo "jemalloc thp:     $(MALLOC_CONF=thp:always,metadata_thp:always LD_PRELOAD=$JE "$PROBE/thp" | grep AnonHuge)"
    echo "jemalloc thp+retain: $(MALLOC_CONF=thp:always,metadata_thp:always,retain:true,dirty_decay_ms:-1,muzzy_decay_ms:-1 LD_PRELOAD=$JE "$PROBE/thp" | grep AnonHuge)"
    if [[ -f $MI ]]; then echo "mimalloc large:   $(MIMALLOC_ALLOW_LARGE_OS_PAGES=1 MIMALLOC_ARENA_EAGER_COMMIT=1 LD_PRELOAD=$MI "$PROBE/thp" | grep AnonHuge)"; fi
    echo "tbbmalloc:        $(LD_LIBRARY_PATH=${TBB%/*}:${LD_LIBRARY_PATH:-} LD_PRELOAD=$TBB "$PROBE/thp" | grep AnonHuge)"
} | tee thp-probe.txt

# --- solver environment (team runtime settings) --------------------------------------
set +eu
module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
module load intel-tbb/2022.3.1   # only puts its lib dir on LD_LIBRARY_PATH; nothing loads it unless preloaded
export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
    WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
set -eu
[[ $WM_OPTIONS == linux64IcxDPInt32OptSPR ]]
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FOAM_SIGFPE=true
export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx I_MPI_DEBUG=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_EXCLUDE_LIST I_MPI_PIN_DOMAIN_SIZE
export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
export I_MPI_PIN_PROCESSOR_LIST
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-1n-104ppn-scatter.dat

variant_env() {   # prints -genv arguments for mpirun (applied to the ranks only)
    case $1 in
        glibc)    ;;
        glibcT)   echo "-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456" ;;
        jemalloc) echo "-genv LD_PRELOAD $JE -genv MALLOC_CONF thp:always,metadata_thp:always" ;;
        # Huge pages without purge/refault churn: keep freed memory (job 180771807: plain
        # thp:always was +45% from re-faulting 2 MiB pages every iteration).
        jemallocR) echo "-genv LD_PRELOAD $JE -genv MALLOC_CONF thp:always,metadata_thp:always,retain:true,dirty_decay_ms:-1,muzzy_decay_ms:-1" ;;
        # Tuning: oversize_threshold:0 keeps >8 MiB allocations (1-node vector fields ~15 MiB/rank)
        # out of jemalloc's huge arena, which purges immediately regardless of dirty_decay_ms;
        # narenas:1 suits single-threaded MPI ranks; tcache_max raises the thread-cache size limit.
        jemallocRO)   echo "-genv LD_PRELOAD $JE -genv MALLOC_CONF thp:always,metadata_thp:always,retain:true,dirty_decay_ms:-1,muzzy_decay_ms:-1,oversize_threshold:0" ;;
        jemallocRON)  echo "-genv LD_PRELOAD $JE -genv MALLOC_CONF thp:always,metadata_thp:always,retain:true,dirty_decay_ms:-1,muzzy_decay_ms:-1,oversize_threshold:0,narenas:1" ;;
        jemallocRONT) echo "-genv LD_PRELOAD $JE -genv MALLOC_CONF thp:always,metadata_thp:always,retain:true,dirty_decay_ms:-1,muzzy_decay_ms:-1,oversize_threshold:0,narenas:1,tcache_max:1048576" ;;
        mimalloc) echo "-genv LD_PRELOAD $MI -genv MIMALLOC_ALLOW_LARGE_OS_PAGES 1 -genv MIMALLOC_ARENA_EAGER_COMMIT 1" ;;
        tbbmalloc) echo "-genv LD_PRELOAD $TBB" ;;
    esac
}

sample_thp() {   # $1 = output file; sums AnonHugePages / Rss over this job's simpleFoam ranks
    sleep 120
    local hp=0 rss=0 n=0 tbb=0 v
    for p in $(pgrep -u "$USER" -x simpleFoam); do
        v=$(awk '/^AnonHugePages/ {h=$2} /^Rss:/ {r=$2} END {print h+0, r+0}' /proc/$p/smaps_rollup 2>/dev/null) || continue
        hp=$((hp + ${v% *})); rss=$((rss + ${v#* })); n=$((n + 1))
        if grep -q '/tbb/2022[.]3/lib/libtbbmalloc[.]so' /proc/$p/maps 2>/dev/null; then tbb=$((tbb + 1)); fi
    done
    echo "ranks=$n AnonHugePages_GiB=$((hp / 1048576)) Rss_GiB=$((rss / 1048576)) ranks_with_tbbmalloc_2022.3=$tbb" > "$1"
}

run() {   # $1 variant
    local v=$1 dir=$ROOT/run_$1
    mkdir "$dir"
    cd "$CASE"
    cp "$TEAM/dicts/1node/controlDict" system/controlDict
    echo "RUN_START variant=$v utc=$(date -u +%FT%TZ)"
    sample_thp "$dir/thp-sample.txt" &
    local sampler=$! status=0
    # shellcheck disable=SC2046
    timeout --kill-after=30s 1200s mpirun -np $NP $(variant_env "$v") simpleFoam -parallel $fh > "$dir/simpleFoam.log" 2>&1 || status=$?
    wait $sampler || true
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    if [[ $status != 0 ]] || ! grep -qx End "$dir/simpleFoam.log" || ! grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log"; then
        echo "RUN_FAILED variant=$v status=$status: $(grep -m1 -iE 'FATAL|error|signal|Segmentation' "$dir/simpleFoam.log" | cut -c1-160)"
        return 0
    fi
    if [[ $v == tbbmalloc ]] && ! grep -q 'ranks_with_tbbmalloc_2022.3=104' "$dir/thp-sample.txt"; then
        echo "RUN_INVALID variant=$v: tbbmalloc 2022.3 not loaded in all ranks ($(cat "$dir/thp-sample.txt"))"
    fi
    local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
    echo "RUN_SUCCESS variant=$v Cd_step40=$cd $(cat "$dir/thp-sample.txt") $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}

# Variant list can be overridden at submission: qsub -v ALLOC_VARIANTS="glibc glibcT jemallocR"
ALLOC_VARIANTS=${ALLOC_VARIANTS:-glibc glibcT jemalloc mimalloc}
echo "VARIANTS $ALLOC_VARIANTS"
for v in $ALLOC_VARIANTS; do run "$v"; done
cp "$TEAM/dicts/1node/controlDict" "$CASE/system/controlDict"

python3 - "$ROOT" $ALLOC_VARIANTS <<'EOF'
import re, sys
from pathlib import Path
root = Path(sys.argv[1])
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
rows, base = [], None
for v in sys.argv[2:]:
    d = root / f'run_{v}'
    found = pat.findall((d / 'simpleFoam.log').read_text()) if (d / 'simpleFoam.log').exists() else []
    thp = (d / 'thp-sample.txt').read_text().strip() if (d / 'thp-sample.txt').exists() else 'n/a'
    t = float(found[-1]) if found else None
    if v == sys.argv[2]: base = t
    rows.append((v, t, thp))
lines = ['# Allocator / huge pages, occDrivAerStaticMesh 65M, 1 node, 104 ranks, 40 steps', '',
         f'ICX OptSPR + Intel MPI, scatter pinning, tuned Allreduce. Reference: {sys.argv[2]} (first variant).', '',
         f'| Variant | s/step | vs {sys.argv[2]} | huge pages held mid-run |', '|---|---:|---:|---|']
for v, t, thp in rows:
    rel = '%+.2f%%' % (100 * (t / base - 1)) if (t and base) else 'n/a'
    lines.append(f"| {v} | {t:.4f} | {rel} | {thp} |" if t else f"| {v} | failed | | {thp} |")
lines += ['', 'Probe on this node:', '```', (root / 'thp-probe.txt').read_text().strip(), '```']
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "ALLOC_COMPLETE root=$ROOT"
