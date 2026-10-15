#!/bin/bash -l
#PBS -N vampir_openfoam
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=32GB
#PBS -l walltime=00:20:00
#PBS -l wd
#PBS -j oe
#PBS -o /home/565/kn3896/scratch/openfoam/logs/profiles/profiling_openfoam_vampir.out

set -eEo pipefail
ROOT=/home/565/kn3896/scratch/openfoam
OPENFOAM_DIR="$ROOT/OpenFOAM-dev"
BASE_CASE="$OPENFOAM_DIR/tutorials/incompressible/icoFoam/cavity/cavity"
PETSC_DIR="$ROOT/ThirdParty-v2606/platforms/linux64IcxDPInt32/petsc-3.25.1"
NPROCS="${PBS_NCPUS:-32}"
MESH_N="${MESH_N:-256}"
RUN_DIR="$ROOT/profiles/vampir_openfoam_${PBS_JOBID:-manual}"
CASE_DIR="$RUN_DIR/case"
PRELOAD_DIR="$RUN_DIR/scorep-preload"
TRACE_DIR="$RUN_DIR/scorep-trace"

mkdir -p "$RUN_DIR" "$PRELOAD_DIR" "$ROOT/logs/profiles"
exec > >(tee "$ROOT/logs/profiles/vampir_openfoam_${PBS_JOBID:-manual}.out") 2>&1
trap 'rc=$?; echo "ERROR: line $LINENO: $BASH_COMMAND (exit $rc)" >&2' ERR

module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2
module load metis/5.1.0 fftw3/3.3.10 intel-mkl/2025.3.1 boost/1.84.0
module load scorep/9.4

export WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx
export WM_MPLIB=INTELMPI I_MPI_CC=icx I_MPI_CXX=icpx
source "$OPENFOAM_DIR/etc/bashrc" || true
[[ "${WM_PROJECT_DIR:-}" == "$OPENFOAM_DIR" ]] || { echo "ERROR: OpenFOAM environment failed"; exit 1; }
[[ "$WM_OPTIONS" == linux64IcxDPInt32Opt ]] || { echo "ERROR: unexpected WM_OPTIONS: $WM_OPTIONS"; exit 1; }

export FOAM_MPI=intelmpi MPI_ARCH_PATH="$I_MPI_ROOT"
export PETSC_ARCH_PATH="$PETSC_DIR" PETSC_ARCH=""
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_THREADING_LAYER=SEQUENTIAL
export LD_LIBRARY_PATH="/apps/scorep/9.4/lib/mpich:$PETSC_DIR/lib:$MKLROOT/lib:${LD_LIBRARY_PATH:-}"
export LIBRARY_PATH="/apps/scorep/9.4/lib/mpich:${LIBRARY_PATH:-}"

cp -a "$BASE_CASE" "$CASE_DIR"
sed -i "s/(20 20 1)/("$MESH_N" "$MESH_N" 1)/" "$CASE_DIR/system/blockMeshDict"
foamDictionary -entry numberOfSubdomains -set "$NPROCS" "$CASE_DIR/system/decomposeParDict"
foamDictionary -entry method -set scotch "$CASE_DIR/system/decomposeParDict"
foamDictionary -entry coeffs -remove "$CASE_DIR/system/decomposeParDict"
foamDictionary -entry deltaT -set 0.0001 "$CASE_DIR/system/controlDict"
foamDictionary -entry endTime -set 0.005 "$CASE_DIR/system/controlDict"
foamDictionary -entry writeInterval -set 50 "$CASE_DIR/system/controlDict"
foamDictionary -entry solvers/p/preconditioner -set distributedDIC "$CASE_DIR/system/fvSolution"
foamDictionary -entry solvers/pFinal/preconditioner -set distributedDIC "$CASE_DIR/system/fvSolution"
foamDictionary -entry solvers/U/solver -set PBiCGStab "$CASE_DIR/system/fvSolution"
foamDictionary -entry solvers/U/smoother -remove "$CASE_DIR/system/fvSolution"
foamDictionary -entry solvers/U/preconditioner -set distributedDILU "$CASE_DIR/system/fvSolution"

blockMesh -case "$CASE_DIR" > "$CASE_DIR/log.blockMesh" 2>&1
decomposePar -case "$CASE_DIR" > "$CASE_DIR/log.decomposePar" 2>&1
PRELOAD_VALUE="$(scorep-preload-init --value-only --mpp=mpi --thread=none --io=posix "$(command -v icoFoam)" "$PRELOAD_DIR")"
[[ -n "$PRELOAD_VALUE" ]] || { echo "ERROR: no Score-P preload value"; exit 1; }
MPI_EVENT_LIB=/apps/scorep/9.4/lib/mpich/libscorep_adapter_mpi_event.so
[[ -r "$MPI_EVENT_LIB" ]] || { echo "ERROR: missing Score-P MPI event adapter: $MPI_EVENT_LIB"; exit 1; }
PRELOAD_VALUE="${PRELOAD_VALUE//-lscorep_adapter_mpi_event/$MPI_EVENT_LIB}"
[[ "$PRELOAD_VALUE" != *-lscorep_adapter_mpi_event* ]] || {
    echo "ERROR: invalid linker option remains in LD_PRELOAD: $PRELOAD_VALUE"
    exit 1
}
echo "Score-P LD_PRELOAD: $PRELOAD_VALUE"

echo "Collecting Score-P OTF2 trace for Vampir on $NPROCS ranks"
/usr/bin/time -v -o "$RUN_DIR/time.txt" mpiexec -n "$NPROCS" env \
    LD_PRELOAD="$PRELOAD_VALUE" LD_LIBRARY_PATH="$LD_LIBRARY_PATH" \
    SCOREP_ENABLE_PROFILING=false SCOREP_ENABLE_TRACING=true \
    SCOREP_EXPERIMENT_DIRECTORY="$TRACE_DIR" SCOREP_TOTAL_MEMORY=1G \
    SCOREP_TRACING_MAX_PROCS="$NPROCS" SCOREP_MPI_ENABLE_GROUPS=CG,ENV \
    icoFoam -parallel -case "$CASE_DIR" > "$CASE_DIR/log.icoFoam" 2>&1

grep -q '^End$' "$CASE_DIR/log.icoFoam" || { echo "ERROR: icoFoam did not finish"; exit 1; }
TRACE_FILE="$TRACE_DIR/traces.otf2"
[[ -f "$TRACE_FILE" ]] || { echo "ERROR: missing OTF2 trace: $TRACE_FILE"; exit 1; }
echo "Vampir trace: $TRACE_FILE"
echo "Trace size: $(du -sh "$TRACE_DIR" | awk '{print $1}')"
