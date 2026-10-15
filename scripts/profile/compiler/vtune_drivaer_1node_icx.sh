#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=01:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N vtune_drivaer_1n
#PBS -j oe

# VTune profile of simpleFoam (custom v2512 ICX OptSPR build + Intel MPI 2021.17.2)
# on occDrivAerStaticMesh 65M, using the team's 1-node setup:
#   /scratch/ya41/apac-2026/openfoam/run/1node/run.sh  (prep chain, hostCollated)
#   /scratch/ya41/apac-2026/openfoam/dicts/1node/      (scotch, 104 ranks, 40 steps)
# plus the team's scatter pinning and tuned Allreduce file.
# Runs: 1) baseline, no tool  2) APS on all ranks  3) VTune on ranks 0 26 52 78.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
NP=104
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_${PBS_JOBID%%.*}
VTUNE_RANKS="0 26 52 78"
fh="-fileHandler hostCollated"

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "VTUNE_JOB_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
cat /proc/sys/kernel/perf_event_paranoid > perf_event_paranoid.txt

# Optional: qsub -v REUSE_CASE=<ROOT of an earlier job> reuses its prepared case (simpleFoam
# runs use noWrite, so the prepared state is not modified) and skips the prep chain.
if [[ -n ${REUSE_CASE:-} ]]; then
    for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$REUSE_CASE/logFiles/$l.log"; done
    ln -s "$REUSE_CASE/case" case
    mkdir logFiles
    cp "$REUSE_CASE"/logFiles/[1-4]*.log logFiles/
    cp "$TEAM/run/1node/intel-scatter-cpus.awk" "$TEAM/run/1node/intelmpi-allreduce-1n-104ppn-scatter.dat" .
    echo "REUSED prepared case $REUSE_CASE/case"
else
# --- case: team case files + linked 65M mesh + team 1-node dicts ----------------
mkdir case logFiles
cp -a "$CASE_SRC/0.orig" "$CASE_SRC/system" case/
mkdir case/constant
cp -a "$CASE_SRC"/constant/{transportProperties,turbulenceProperties} case/constant/
ln -s "$MESH" case/constant/polyMesh
cp "$TEAM/dicts/1node/caseDefinition"   case/system/include/caseDefinition
# dicts/1node/decomposeParDict references $nDomains/$nDomainWeights, which dicts/1node/caseDefinition
# never defines (decomposePar aborts). decomposeParDict_base is the same file without that
# scotchCoeffs block: plain scotch into $nCores=104 parts.
cp "$TEAM/dicts/decomposeParDict_base" case/system/decomposeParDict
cp "$TEAM/run/1node/intel-scatter-cpus.awk" "$TEAM/run/1node/intelmpi-allreduce-1n-104ppn-scatter.dat" .

# --- prep chain exactly as run.sh (site module) ---------------------------------
(
    set +eu
    module purge
    module load OpenFOAM/2512
    . "$WM_PROJECT_DIR/bin/tools/RunFunctions"
    set -eu
    cd case
    par=(mpirun -np $NP)
    decomposePar -constant $fh > ../logFiles/10_decomposePar.log 2>&1
    restore0Dir -processor > ../logFiles/15_restore0Dir.log 2>&1
    "${par[@]}" renumberMesh -constant -overwrite -parallel $fh > ../logFiles/20_renumberMesh.log 2>&1
    "${par[@]}" potentialFoam -initialiseUBCs -parallel $fh > ../logFiles/30_potentialFoam.log 2>&1
    "${par[@]}" applyBoundaryLayer -ybl "0.0450244" -parallel $fh > ../logFiles/40_applyBoundaryLayer.log 2>&1
    for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' ../logFiles/$l.log; done
)
echo "PREPARED $(grep -m1 'Number of cells' logFiles/10_decomposePar.log || true)"
fi

# dicts/1node/controlDict is a full controlDict.noWrite with endTime 40 / runTimeModifiable no.
# run.sh copies it to system/include/controlDict, which system/controlDict never includes, so
# simpleFoam would run the upstream endTime 4000. Install it as system/controlDict instead.
cp "$TEAM/dicts/1node/controlDict" case/system/controlDict
[[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' case/system/controlDict) == 40 ]]

# --- ICX OptSPR + Intel MPI environment (team runtime settings) -----------------
set +eu
module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
export FOAM_EXTRA_CFLAGS="-xSAPPHIRERAPIDS -fp-model=precise" FOAM_EXTRA_CXXFLAGS="-xSAPPHIRERAPIDS -fp-model=precise"
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
I_MPI_PIN_PROCESSOR_LIST=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f ./intel-scatter-cpus.awk)
export I_MPI_PIN_PROCESSOR_LIST
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-1n-104ppn-scatter.dat
{
    echo "simpleFoam=$(command -v simpleFoam)"
    echo "compiler=$(icpx --version | head -1)"
    echo "mpi=$(mpirun --version 2>&1 | head -1)"
    echo "vtune=$(vtune --version | head -1)"
    echo "pin_list=$I_MPI_PIN_PROCESSOR_LIST"
} > build-metadata.txt

# Pick the richest VTune analysis this node allows (user-mode only when paranoid=2).
ANALYSIS=
# Job 180728663: driverless hpc-performance hung (amplxe-perf never exited after a 2 s sleep),
# so every self-test is time-limited and user-mode (software) sampling is the perf-free fallback.
ANALYSIS_OPTS=
for cfg in "hpc-performance:" "hotspots:-knob sampling-mode=hw" "hotspots:-knob sampling-mode=sw"; do
    a=${cfg%%:*} opts=${cfg#*:} tag=$a-${opts##*=}
    if timeout --kill-after=10s 120s vtune -collect $a $opts -r "$ROOT/vtune-selftest-$tag" -- /bin/sleep 2 \
            > "vtune-selftest-$tag.log" 2>&1; then
        ANALYSIS=$a ANALYSIS_OPTS=$opts; break
    fi
    echo "VTune self-test failed or timed out: $a $opts"
    pkill -u "$USER" -f amplxe-perf || true
done
[[ -n $ANALYSIS ]] || { echo "ERROR: no VTune analysis works here (see vtune-selftest-*.log)"; exit 1; }
echo "VTUNE_ANALYSIS=$ANALYSIS $ANALYSIS_OPTS"

cat > vtune-wrap.sh <<EOF
#!/bin/bash
case " $VTUNE_RANKS " in
    *" \$PMI_RANK "*) exec vtune -collect $ANALYSIS $ANALYSIS_OPTS -data-limit=0 -r "$ROOT/vtune_rank\$PMI_RANK" -- "\$@" ;;
    *) exec "\$@" ;;
esac
EOF
chmod +x vtune-wrap.sh

run() {   # $1 label, rest = launcher prefix after mpirun -np
    local label=$1; shift
    echo "RUN_START $label utc=$(date -u +%FT%TZ)"
    (cd case && timeout --kill-after=30s 1500s mpirun -np $NP "$@" simpleFoam -parallel $fh > "$ROOT/logFiles/50_simpleFoam_$label.log" 2>&1)
    grep -qx End "logFiles/50_simpleFoam_$label.log"
    echo "RUN_SUCCESS $label $(grep 'Average wall-clock time per time step' logFiles/50_simpleFoam_$label.log | tail -1)"
}

run baseline
run aps aps --result-dir="$ROOT/aps_result"
run vtune "$ROOT/vtune-wrap.sh"

# --- reports ------------------------------------------------------------------
aps --report="$ROOT/aps_result" > aps-report.txt 2>&1 || true
for r in $VTUNE_RANKS; do
    vtune -report summary -r "vtune_rank$r" > "vtune-summary-rank$r.txt" 2>&1 || true
    vtune -report hotspots -r "vtune_rank$r" -group-by function -format csv -csv-delimiter comma \
        > "vtune-hotspots-rank$r.csv" 2>&1 || true
done
grep -h 'Average wall-clock time per time step' logFiles/50_simpleFoam_*.log
python3 - "$ROOT" <<'EOF'
import csv, sys
from pathlib import Path
root = Path(sys.argv[1])
rows = list(csv.reader((root / 'vtune-hotspots-rank0.csv').open()))
if rows:
    print('Top functions, rank 0:')
    for r in rows[1:16]:
        print('  ' + ' | '.join(r[:3]))
EOF
echo "VTUNE_JOB_COMPLETE root=$ROOT analysis=$ANALYSIS $ANALYSIS_OPTS"
