#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N vtune_cur_1n
#PBS -j oe

# VTune hotspots of the current Bonus setup (measured 4n: simpleFoamFused + fixedIterSolvers v2 +
# fusedVectorSolve v3 + stepReuse, cache grad(U), grad(p), grad(k), grad(omega)), to rank what is left.
# 1 node x 104 ranks, prepared 1-node case (no prep), 20 steps. VTune user-mode (software) sampling on
# ranks 0 26 52 78; hotspots, top-down (inclusive) and memcpy/memset callers.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
CASE=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/vtune/current_1node_${PBS_JOBID%%.*}
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
RLIB=/scratch/ya41/kn3896/openfoam/tools/stepReuse/lib/linux64IcxDPInt32OptSPR/libstepReuse.so
SFF=/scratch/ya41/kn3896/openfoam/tools/simpleFoamFused/bin/linux64IcxDPInt32OptSPR/simpleFoamFused
VTUNE_RANKS="0 26 52 78"
STEPS=20
fh=(-fileHandler hostCollated)

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "VTUNE_CURRENT_FAILED line=$LINENO status=$?"; cp "$ROOT/fvSolution.orig" "$CASE/system/fvSolution" 2>/dev/null || true' ERR
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
    *" \$PMI_RANK "*) exec vtune -collect hotspots -knob sampling-mode=sw -data-limit=0 -r "$ROOT/vtune_rank\$PMI_RANK" -- env LD_PRELOAD="$FLIB:$RLIB" "\$@" ;;
    *) exec env LD_PRELOAD="$FLIB:$RLIB" "\$@" ;;
esac
EOF
chmod +x vtune-wrap.sh

cd "$CASE"
cp "$TEAM/dicts/1node/controlDict" system/controlDict
sed -i "s/^endTime .*/endTime           $STEPS;/" system/controlDict
cp "$ROOT/fvSolution.fixed" system/fvSolution
echo "libs (\"$LIB\");" >> system/controlDict
echo "RUN_START utc=$(date -u +%FT%TZ)"
status=0
timeout --kill-after=30s 1500s mpirun -np 104 "${glibcT[@]}" -genv FUSED_VECTOR_SOLVE_REPORT 1 -genv STEP_REUSE_REPORT 1 -genv SIMPLEFOAM_FUSED_REPORT 1 "$ROOT/vtune-wrap.sh" "$SFF" -parallel "${fh[@]}" > "$ROOT/simpleFoam.log" 2>&1 || status=$?
if [[ -d postProcessing ]]; then mv postProcessing "$ROOT/"; fi
cp "$ROOT/fvSolution.orig" system/fvSolution
cp "$TEAM/dicts/1node/controlDict" system/controlDict
cd "$ROOT"
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"   # case restored
[[ $status == 0 ]]
grep -qx End simpleFoam.log
echo "RUN_SUCCESS $(grep 'Average wall-clock time per time step' simpleFoam.log | tail -1)"

grep -q 'fusedVectorSolve: U fused path' simpleFoam.log
grep -q 'simpleFoamFused: pressure equation steps fused' simpleFoam.log
grep -q 'stepReuse: laplacian: fused' simpleFoam.log
for r in $VTUNE_RANKS; do
    d=$(ls -d vtune_rank$r.* | head -1)   # VTune appends the host name
    vtune -report summary -r "$d" > "vtune-summary-rank$r.txt" 2>&1 || true
    vtune -report hotspots -r "$d" -group-by function -format csv -csv-delimiter comma \
        > "vtune-hotspots-rank$r.csv" 2>&1 || true
    vtune -report top-down -r "$d" -format csv -csv-delimiter comma \
        > "vtune-topdown-rank$r.csv" 2>&1 || true
done
for f in __memcpy_avx512_unaligned_erms __intel_avx_rep_memset __memset_avx512_unaligned_erms; do
    vtune -report callstacks -r "$(ls -d vtune_rank26.* | head -1)" -filter function="$f" -format csv -csv-delimiter comma \
        > "callers-$f-rank26.csv" 2>&1 || true
done
python3 - "$ROOT" <<'PY2'
import csv, sys
from pathlib import Path
root = Path(sys.argv[1])
rows = [r for r in csv.reader((root / 'vtune-hotspots-rank26.csv').open()) if r][1:]
tot = sum(float(r[1]) for r in rows if r[1].replace('.', '', 1).isdigit())
print(f'Top functions, rank 26 (total CPU {tot:.1f} s):')
for r in rows[:20]:
    print(f'  {100 * float(r[1]) / tot:5.1f}%  {r[0][:90]}')
PY2
echo "VTUNE_CURRENT_COMPLETE root=$ROOT"
