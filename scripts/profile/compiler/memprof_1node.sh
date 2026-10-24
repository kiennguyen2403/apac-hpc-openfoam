#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb,jobfs=20gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N memprof_1n
#PBS -j oe

# Memory-stall profile of the current Bonus setup (simpleFoamFused + fixedIterSolvers v2 + fusedVectorSolve v3 +
# stepReuse (climit), cache grad(U/p/k/omega), UCX zero-copy off), 1 node x 104 ranks, prepared 1-node case.
#   aps   - APS on all ranks, 6 steps: memory/DRAM stall shares (October, original setup: 67.7% memory stalls,
#           49.3% DRAM stalls; job 180732331)
#   uarch - VTune uarch-exploration (driverless perf, user mode) on ranks 26 and 78 while the other ranks run
#           normally, 6 steps: per function Memory Bound, L1/L2/L3/DRAM Bound, DRAM Bandwidth vs Latency.
# Goal: is what remains DRAM-bandwidth bound (fix: fewer bytes) or latency bound (fix: huge pages, prefetch)?

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
CASE=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/vtune/memprof_1node_${PBS_JOBID%%.*}
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
RLIB=/scratch/ya41/kn3896/openfoam/tools/stepReuse/lib/linux64IcxDPInt32OptSPR/libstepReuse.so
SFF=/scratch/ya41/kn3896/openfoam/tools/simpleFoamFused/bin/linux64IcxDPInt32OptSPR/simpleFoamFused
VTUNE_RANKS="26 78"
STEPS=6
fh=(-fileHandler hostCollated)

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "MEMPROF_FAILED line=$LINENO status=$?"; cp "$ROOT/fvSolution.orig" "$CASE/system/fvSolution" 2>/dev/null || true' ERR
cd "$ROOT"
cp "$TEAM/run/1node/intel-scatter-cpus.awk" "$TEAM/run/1node/intelmpi-allreduce-1n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/../logFiles/$l.log"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
ls "$CASE/dynamicCode/platforms/linux64IcxDPInt32OptSPR/lib/"libwallClockTimeStatistics_*.so > /dev/null
[[ -f $LIB && -f $FLIB && -f $RLIB && -x $SFF ]]
for l in "$LIB" "$FLIB" "$RLIB" "$SFF"; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
cp "$CASE/system/fvSolution" fvSolution.orig
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.orig > fvSolution.fixed
echo "fvSolution changes (base -> fixed):"; diff fvSolution.orig fvSolution.fixed || true
[[ $(grep -c 'smoothSolverFixed;' fvSolution.fixed) == 3 && $(grep -c 'PCGFixed;' fvSolution.fixed) == 1 ]]
awk '{print} /^[[:space:]]*grad\(U\);/ && incache {print "    grad(p);"; print "    grad(k);"; print "    grad(omega);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.fixed > fvSolution.tmp && mv fvSolution.tmp fvSolution.fixed
[[ $(grep -c 'grad(p);' fvSolution.fixed) == 1 && $(grep -c 'grad(omega);' fvSolution.fixed) == 1 ]]
if [[ -d $CASE/postProcessing ]]; then mv "$CASE/postProcessing" "$ROOT/postProcessing.before"; fi

set +eu
module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
    WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
module load intel-vtune/2025.8.1
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
glibcT=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)

cat > vtune-wrap.sh <<EOF
#!/bin/bash
case " $VTUNE_RANKS " in
    *" \$PMI_RANK "*) exec vtune -collect uarch-exploration -data-limit=0 -r "$ROOT/uarch_rank\$PMI_RANK" -- env LD_PRELOAD="$FLIB:$RLIB" "\$@" ;;
    *) exec env LD_PRELOAD="$FLIB:$RLIB" "\$@" ;;
esac
EOF
chmod +x vtune-wrap.sh

cd "$CASE"
cp "$TEAM/dicts/1node/controlDict" system/controlDict
sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
cp "$ROOT/fvSolution.fixed" system/fvSolution
echo "libs (\"$LIB\");" >> system/controlDict
ucxT=(-genv UCX_ZCOPY_THRESH inf -genv UCX_RNDV_THRESH inf)
envs=(-genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1 -genv SIMPLEFOAM_FUSED_REPORT 1)

echo "RUN_START aps utc=$(date -u +%FT%TZ)"
status=0
timeout --kill-after=30s 1200s mpirun -np 104 "${glibcT[@]}" "${ucxT[@]}" -genv LD_PRELOAD "$FLIB:$RLIB" "${envs[@]}" \
    aps --result-dir="$ROOT/aps_result" "$SFF" -parallel "${fh[@]}" > "$ROOT/simpleFoam_aps.log" 2>&1 || status=$?
grep -qx End "$ROOT/simpleFoam_aps.log" || { echo "APS run failed status=$status"; false; }
[[ $status == 0 ]] || echo "WARNING aps: mpirun status $status after End (APS shutdown); data kept"
echo "RUN_SUCCESS aps $(grep 'Average wall-clock time per time step' "$ROOT/simpleFoam_aps.log" | tail -1)"
if [[ -d postProcessing ]]; then rm -rf postProcessing; fi

echo "RUN_START uarch utc=$(date -u +%FT%TZ)"
status=0
timeout --kill-after=30s 1500s mpirun -np 104 "${glibcT[@]}" "${ucxT[@]}" "${envs[@]}" \
    "$ROOT/vtune-wrap.sh" "$SFF" -parallel "${fh[@]}" > "$ROOT/simpleFoam_uarch.log" 2>&1 || status=$?
if [[ -d postProcessing ]]; then rm -rf postProcessing; fi
cp "$ROOT/fvSolution.orig" system/fvSolution
cp "$TEAM/dicts/1node/controlDict" system/controlDict
cd "$ROOT"
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"   # case restored
[[ $status == 0 ]]
grep -qx End simpleFoam_uarch.log
echo "RUN_SUCCESS uarch $(grep 'Average wall-clock time per time step' simpleFoam_uarch.log | tail -1)"
grep -q 'simpleFoamFused: pressure equation steps fused' simpleFoam_uarch.log

aps-report "$ROOT/aps_result" > aps-summary.txt 2>&1 || true
aps-report -m "$ROOT/aps_result" > aps-metrics.txt 2>&1 || true
for r in $VTUNE_RANKS; do
    d=$(ls -d uarch_rank$r.* | head -1)   # VTune appends the host name
    vtune -report summary -r "$d" > "uarch-summary-rank$r.txt" 2>&1 || true
    vtune -report hotspots -r "$d" -group-by function -format csv -csv-delimiter comma \
        > "uarch-functions-rank$r.csv" 2>&1 || true
done
grep -E "Memory Stalls|DRAM Stalls|Cache Stalls|Memory Bound|DRAM Bound|Bandwidth|Latency" aps-summary.txt | head -12 || true
grep -E "Memory Bound|DRAM Bound|Memory Bandwidth|Memory Latency|L1 Bound|L2 Bound|L3 Bound|Store Bound|Elapsed|CPI" uarch-summary-rank26.txt | head -20 || true
echo "MEMPROF_COMPLETE root=$ROOT"
