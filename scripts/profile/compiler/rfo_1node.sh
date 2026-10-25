#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=00:25:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N rfo_1n
#PBS -j oe

# Store-side memory traffic of the current Bonus setup (with compactAddr), 1 node x 104 ranks, prepared 1-node
# case, 12 steps: how much of the L2 fill traffic is read-for-ownership (the fills streaming stores could
# avoid). perf stat -I 1000 user mode on every rank: L2_LINES_IN.ALL, L2_RQSTS.ALL_RFO / RFO_MISS / HWPF_MISS,
# OFFCORE_REQUESTS.DEMAND_RFO, L2_LINES_OUT.NON_SILENT (dirty write-backs = upper bound of write-only fills).
# tools/memBandwidth/analyse_rfo.py on the time-step window. Memory roofline of the same setup: job 182171696.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
CASE=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
RLIB=/scratch/ya41/kn3896/openfoam/tools/stepReuse/lib/linux64IcxDPInt32OptSPR/libstepReuse.so
SFF=/scratch/ya41/kn3896/openfoam/tools/simpleFoamFused/bin/linux64IcxDPInt32OptSPR/simpleFoamFused
MBW=/scratch/ya41/kn3896/openfoam/tools/memBandwidth
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/vtune/rfo_1node_${PBS_JOBID%%.*}
CALIB=/scratch/ya41/kn3896/openfoam/tools/compactAddr/lib/linux64IcxDPInt32OptSPR/libcompactAddr.so
STEPS=12
fh=(-fileHandler hostCollated)

mkdir -p "$ROOT/perf"
exec > >(tee "$ROOT/job.log") 2>&1
restore() {
    cp "$ROOT/fvSolution.orig" "$CASE/system/fvSolution" 2>/dev/null || true
    cp "$TEAM/dicts/1node/controlDict" "$CASE/system/controlDict" 2>/dev/null || true
    if [[ -d $CASE/postProcessing ]]; then mv "$CASE/postProcessing" "$ROOT/postProcessing.run"; fi
}
trap 'echo "RFO_FAILED line=$LINENO status=$?"; restore' ERR
cd "$ROOT"
cp "$TEAM/run/1node/intel-scatter-cpus.awk" "$TEAM/run/1node/intelmpi-allreduce-1n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/../logFiles/$l.log"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
[[ ! -d $CASE/postProcessing ]]
[[ -f $LIB && -f $FLIB && -f $RLIB && -f $CALIB && -x $SFF ]]
for l in "$LIB" "$FLIB" "$RLIB" "$CALIB" "$SFF"; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
cp "$CASE/system/fvSolution" fvSolution.orig
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.orig > fvSolution.fixed
[[ $(grep -c 'smoothSolverFixed;' fvSolution.fixed) == 3 && $(grep -c 'PCGFixed;' fvSolution.fixed) == 1 ]]
awk '{print} /^[[:space:]]*grad\(U\);/ && incache {print "    grad(p);"; print "    grad(k);"; print "    grad(omega);"} /^cache/ {incache=1} /^}/ {incache=0}' fvSolution.fixed > fvSolution.tmp && mv fvSolution.tmp fvSolution.fixed
[[ $(grep -c 'grad(p);' fvSolution.fixed) == 1 && $(grep -c 'grad(omega);' fvSolution.fixed) == 1 ]]

set +eu
module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
    WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
set -eu
[[ $WM_OPTIONS == linux64IcxDPInt32OptSPR ]]
command -v perf
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FOAM_SIGFPE=true
export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx I_MPI_DEBUG=0
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_EXCLUDE_LIST I_MPI_PIN_DOMAIN_SIZE
export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
export I_MPI_PIN_PROCESSOR_LIST
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-1n-104ppn-scatter.dat
glibcT=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)

# 2) solver with per-rank perf counters (user mode)
EV=cycles:u,instructions:u
EV+=,cpu/event=0x25,umask=0x1f,name=l2_lines_in/u
EV+=,cpu/event=0x24,umask=0xe2,name=rfo_all/u
EV+=,cpu/event=0x24,umask=0x22,name=rfo_miss/u
EV+=,cpu/event=0x21,umask=0x04,name=offcore_rfo/u
EV+=,cpu/event=0x24,umask=0x30,name=hwpf_miss/u
EV+=,cpu/event=0x26,umask=0x02,name=wb_dirty/u
# every event must be accepted, otherwise ranks would fail and the rest hang (as VTune did in 182170650)
# (4 general counters per hyperthread, 8 events: perf multiplexes, so test with a few seconds of busy work and
#  fail only on events the CPU rejects; "<not counted>" in a short test is just multiplexing)
perf stat -x, -e "$EV" -o perf-selftest.txt -- bash -c 'i=0; while ((i < 4000000)); do ((i++)); done'
cat perf-selftest.txt
[[ $(grep -ciE 'not supported|invalid|unknown' perf-selftest.txt || true) == 0 ]]
cat > perf-wrap.sh <<WRAP
#!/bin/bash
exec perf stat -I 1000 -x, -o "$ROOT/perf/rank\$PMI_RANK.txt" -e "$EV" -- env LD_PRELOAD="$CALIB:$FLIB:$RLIB" "\$@"
WRAP
chmod +x perf-wrap.sh

cd "$CASE"
cp "$TEAM/dicts/1node/controlDict" system/controlDict
sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
cp "$ROOT/fvSolution.fixed" system/fvSolution
echo "libs (\"$LIB\");" >> system/controlDict
echo "RUN_START utc=$(date -u +%FT%TZ)"
status=0
timeout --kill-after=30s 900s mpirun -np 104 "${glibcT[@]}" \
    -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1 -genv SIMPLEFOAM_FUSED_REPORT 1 \
    "$ROOT/perf-wrap.sh" "$SFF" -parallel "${fh[@]}" > "$ROOT/simpleFoam.log" 2>&1 || status=$?
cd "$ROOT"
restore
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"   # case restored
[[ -z $(ls -d "$CASE"/processors104/[1-9]* 2>/dev/null) ]]             # nothing written
[[ $status == 0 ]]
grep -qx End simpleFoam.log
grep -q 'simpleFoamFused: pressure equation steps fused' simpleFoam.log
grep -q 'compactAddr: GaussSeidel: 16-bit neighbour offsets' simpleFoam.log
echo "RUN_SUCCESS $(grep 'Average wall-clock time per time step' simpleFoam.log | tail -1)"
[[ $(ls perf/rank*.txt | wc -l) == 104 ]]

# 3) analysis of the time-step window
python3 "$MBW/analyse_rfo.py" "$ROOT/perf" "$ROOT/simpleFoam.log" | tee analysis.txt
n=$(find perf | wc -l); tar czf perf.tar.gz perf && [[ $(tar tzf perf.tar.gz | wc -l) == "$n" ]] && rm -rf perf
echo "RFO_COMPLETE root=$ROOT"
