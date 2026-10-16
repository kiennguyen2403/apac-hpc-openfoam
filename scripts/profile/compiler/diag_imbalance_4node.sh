#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb,jobfs=40gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N diag_imbal_4n
#PBS -j oe

# Load-imbalance diagnosis with the current Bonus setup (fixedIterSolvers plugin + fused vector
# solve v2, glibcT, ICX OptSPR + Intel MPI, scatter pinning, tuned Allreduce), 4 nodes x 104,
# scotch (4 104)/(1 0.5). Per-rank compute = APS elapsed - MPI time.
#   normal5  : normal mapping, 5 steps     } per-step compute per rank = normal25 - normal5
#   normal25 : normal mapping, 25 steps    }
#   perm25   : permuted mapping, 25 steps: node order reversed and the core list reversed inside
#              each node, so every subdomain runs on a different node, socket and core.
#              Slow ranks that stay slow -> subdomain cause; slowness that moves -> placement.
# jobfs=40gb: APS stages its data in the job-local directory (100 MB default per node killed job
# 181350412). I_MPI_DEBUG=4 logs the rank -> node/cpu map of every run. APS reports are generated afterwards
# on the login node (aps-report -t). Runs in a linked copy of case416_w05.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
LIB=/scratch/ya41/kn3896/openfoam/tools/fixedIterSolvers/lib/linux64IcxDPInt32OptSPR/libfixedIterSolvers.so
FLIB=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so
CUR=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/diag/imbalance_4node_${PBS_JOBID%%.*}
CASEDIR=$ROOT/case
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "DIAG_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]]
cp "$TEAM/run/4node/intel-scatter-cpus.awk" "$TEAM/run/4node/intelmpi-allreduce-4n-104ppn-scatter.dat" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CUR/logFiles/$l.log"; done
[[ -f $LIB && -f $FLIB ]]
cp "$CUR/system/fvSolution.fixedIter" fvSolution.basic
sed -e 's/^\(\s*solver\s\+\)smoothSolver;/\1smoothSolverFixed;/' -e 's/^\(\s*solver\s\+\)PCG;/\1PCGFixed;/' fvSolution.basic > fvSolution.bonus
[[ $(grep -c 'smoothSolverFixed;' fvSolution.bonus) == 3 && $(grep -c 'PCGFixed;' fvSolution.bonus) == 1 ]]
lscpu -p=CPU,CORE,SOCKET,NODE > lscpu.txt

mkdir -p "$CASEDIR/constant"
cp -a "$CUR/system" "$CASEDIR/"
cp fvSolution.bonus "$CASEDIR/system/fvSolution"
cp -a "$CUR"/constant/{transportProperties,turbulenceProperties} "$CASEDIR/constant/"
ln -s "$MESH" "$CASEDIR/constant/polyMesh"
for p in "$CUR"/processors416_*; do ln -s "$p" "$CASEDIR/${p##*/}"; done
cp -a "$CUR/dynamicCode" "$CASEDIR/dynamicCode"

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
export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx I_MPI_DEBUG=4
unset I_MPI_PIN_DOMAIN I_MPI_PIN_PROCESSOR_EXCLUDE_LIST I_MPI_PIN_DOMAIN_SIZE
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
export I_MPI_PIN=1 I_MPI_PIN_ORDER=scatter I_MPI_PIN_CELL=core I_MPI_HYDRA_TOPOLIB=hwloc
PIN_NORMAL=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
PIN_PERM=$(tr ',' '\n' <<< "$PIN_NORMAL" | tac | paste -sd,)
[[ $(tr ',' '\n' <<< "$PIN_NORMAL" | wc -l) == 104 && $(tr ',' '\n' <<< "$PIN_PERM" | wc -l) == 104 ]]
HOSTS_NORMAL=$(sort -u "$PBS_NODEFILE" | paste -sd,)
HOSTS_PERM=$(sort -u "$PBS_NODEFILE" | tac | paste -sd,)
echo "hosts normal: $HOSTS_NORMAL"; echo "hosts perm:   $HOSTS_PERM"
unset I_MPI_ADJUST_ALLREDUCE I_MPI_ADJUST_ALLREDUCE_LIST I_MPI_TUNING_MODE
export I_MPI_TUNING_BIN=$ROOT/intelmpi-allreduce-4n-104ppn-scatter.dat
glibcT=(-genv MALLOC_MMAP_THRESHOLD_ 33554432 -genv MALLOC_TRIM_THRESHOLD_ 4294967296 -genv MALLOC_TOP_PAD_ 268435456)

run() {   # $1 name, $2 steps, $3 hosts, $4 pin list
    local name=$1 steps=$2 dir=$ROOT/run_$1 status=0
    mkdir "$dir"
    cd "$CASEDIR"
    cp "$TEAM/dicts/4node/controlDict" system/controlDict
    sed -i "s/^endTime .*/endTime           $steps;/" system/controlDict
    echo "libs (\"$LIB\");" >> system/controlDict
    echo "RUN_START $name steps=$steps utc=$(date -u +%FT%TZ)"
    I_MPI_PIN_PROCESSOR_LIST=$4 timeout --kill-after=30s 900s \
        mpirun -np 416 -ppn 104 -hosts "$3" "${glibcT[@]}" -genv LD_PRELOAD "$FLIB" -genv FUSED_VECTOR_SOLVE_REPORT 1 \
        aps --result-dir="$dir/aps_result" simpleFoam -parallel "${fh[@]}" > "$dir/simpleFoam.log" 2>&1 || status=$?
    if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
    cd "$ROOT"
    # Job 181350279: the run completed and APS wrote its data, but the hydra proxy on one node
    # exited abnormally at shutdown (APS finalisation), giving a non-zero mpirun status.
    if [[ $status != 0 ]]; then
        if grep -qx End "$dir/simpleFoam.log" && [[ -d $dir/aps_result ]]; then
            echo "WARNING $name: mpirun status $status after End (APS shutdown); APS data kept"
        else
            false
        fi
    fi
    grep -qx End "$dir/simpleFoam.log"
    [[ $(grep -c '^Time = ' "$dir/simpleFoam.log") == $steps ]]
    grep -q 'fusedVectorSolve: U fused path' "$dir/simpleFoam.log"
    echo "RUN_SUCCESS $name pinmap_lines=$(grep -c -E '^\[0\] MPI startup\(\): [0-9]+ +[0-9]+ +' "$dir/simpleFoam.log") $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
}
run normal5  5  "$HOSTS_NORMAL" "$PIN_NORMAL"
run normal25 25 "$HOSTS_NORMAL" "$PIN_NORMAL"
run perm25   25 "$HOSTS_PERM"   "$PIN_PERM"
echo "DIAG_COMPLETE root=$ROOT"
