#!/bin/bash -l
#PBS -N test_openfoam
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=8GB
#PBS -l walltime=00:15:00
#PBS -l wd
#PBS -j oe
#PBS -o /home/565/kn3896/scratch/openfoam/logs/tests/test_openfoam.out

set -eEo pipefail

OPENFOAM_ROOT=/home/565/kn3896/scratch/openfoam
OPENFOAM_DIR="$OPENFOAM_ROOT/OpenFOAM-dev"
THIRDPARTY_DIR="$OPENFOAM_ROOT/ThirdParty-v2606"
TUTORIAL="$OPENFOAM_DIR/tutorials/IO/cavity_parProfiling"
PETSC_DIR="$THIRDPARTY_DIR/platforms/linux64IcxDPInt32/petsc-3.25.1"
NPROCS="${PBS_NCPUS:-4}"
RUN_DIR="$OPENFOAM_ROOT/tests/openfoam_intel_mkl_${PBS_JOBID:-manual}"
LOG_DIR="$OPENFOAM_ROOT/logs/tests"

mkdir -p "$RUN_DIR" "$LOG_DIR"
exec > >(tee "$LOG_DIR/test_openfoam_${PBS_JOBID:-manual}.out") 2>&1
trap 'rc=$?; echo "ERROR: line $LINENO: $BASH_COMMAND (exit $rc)" >&2' ERR

module purge
module load intel-compiler-llvm/2025.3.2
module load intel-mpi/2021.17.2
module load metis/5.1.0
module load fftw3/3.3.10
module load intel-mkl/2025.3.1
module load boost/1.84.0

export WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx
export WM_MPLIB=INTELMPI I_MPI_CC=icx I_MPI_CXX=icpx
source "$OPENFOAM_DIR/etc/bashrc" || true

[[ "${WM_PROJECT_DIR:-}" == "$OPENFOAM_DIR" ]] || { echo "ERROR: wrong OpenFOAM environment: ${WM_PROJECT_DIR:-unset}"; exit 1; }
[[ "$WM_OPTIONS" == linux64IcxDPInt32Opt ]] || { echo "ERROR: unexpected WM_OPTIONS: $WM_OPTIONS"; exit 1; }

export FOAM_MPI=intelmpi MPI_ARCH_PATH="$I_MPI_ROOT"
export PETSC_ARCH_PATH="$PETSC_DIR" PETSC_ARCH=""
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_THREADING_LAYER=SEQUENTIAL
export LD_LIBRARY_PATH="$PETSC_DIR/lib:$MKLROOT/lib:${LD_LIBRARY_PATH:-}"

require_command()
{
    command -v "$1" >/dev/null 2>&1 || { echo "ERROR: required command not found: $1"; exit 1; }
}

check_log()
{
    local log="$1"
    grep -q '^End$' "$log" || { echo "ERROR: successful termination not found in $log"; exit 1; }
    ! grep -Eq 'FOAM FATAL|Segmentation fault' "$log" || { echo "ERROR: failure found in $log"; exit 1; }
}

for cmd in blockMesh decomposePar icoFoam foamDictionary mpiexec; do require_command "$cmd"; done
[[ -f "$PETSC_DIR/lib/libpetsc.so" ]] || { echo "ERROR: PETSc library not found"; exit 1; }

PETSC_LDD="$(ldd "$PETSC_DIR/lib/libpetsc.so")"
! grep -q 'not found' <<< "$PETSC_LDD" || { echo "$PETSC_LDD"; exit 1; }
for dependency in libmkl_core libmkl_intel_lp64 libmkl_sequential libHYPRE libmpi; do
    grep -q "$dependency" <<< "$PETSC_LDD" || { echo "ERROR: PETSc is not linked to $dependency"; exit 1; }
done

echo "OpenFOAM: $WM_PROJECT_DIR"
echo "Build: $WM_OPTIONS"
echo "PETSc: $PETSC_DIR"
echo "MPI ranks: $NPROCS"

cp -a "$TUTORIAL" "$RUN_DIR/serial"
cp -a "$TUTORIAL" "$RUN_DIR/parallel"

cp "$RUN_DIR/serial/system/fvSolution.template" "$RUN_DIR/serial/system/fvSolution"
blockMesh -case "$RUN_DIR/serial" > "$RUN_DIR/serial/log.blockMesh" 2>&1
icoFoam -case "$RUN_DIR/serial" > "$RUN_DIR/serial/log.icoFoam" 2>&1
check_log "$RUN_DIR/serial/log.icoFoam"

cp "$RUN_DIR/parallel/system/fvSolution.template" "$RUN_DIR/parallel/system/fvSolution"
foamDictionary -entry numberOfSubdomains -set "$NPROCS" "$RUN_DIR/parallel/system/decomposeParDict"
foamDictionary -entry method -set scotch "$RUN_DIR/parallel/system/decomposeParDict"
foamDictionary -entry coeffs -remove "$RUN_DIR/parallel/system/decomposeParDict"
blockMesh -case "$RUN_DIR/parallel" > "$RUN_DIR/parallel/log.blockMesh" 2>&1
decomposePar -case "$RUN_DIR/parallel" > "$RUN_DIR/parallel/log.decomposePar" 2>&1
mpiexec -n "$NPROCS" icoFoam -parallel -case "$RUN_DIR/parallel" > "$RUN_DIR/parallel/log.icoFoam" 2>&1
check_log "$RUN_DIR/parallel/log.icoFoam"

echo
echo "OpenFOAM Intel/MKL tests passed"
echo "Results: $RUN_DIR"
