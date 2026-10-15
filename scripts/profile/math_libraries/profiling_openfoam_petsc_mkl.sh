#!/bin/bash -l
#PBS -N profile_petsc_mkl
#PBS -q normal
#PBS -l ncpus=8
#PBS -l mem=32GB
#PBS -l walltime=00:30:00
#PBS -l wd
#PBS -j oe
#PBS -o /home/565/kn3896/scratch/openfoam/logs/profiles/profiling_petsc_mkl.out

set -eEo pipefail

OPENFOAM_ROOT=/home/565/kn3896/scratch/openfoam
OPENFOAM_DIR="$OPENFOAM_ROOT/OpenFOAM-dev"
THIRDPARTY_DIR="$OPENFOAM_ROOT/ThirdParty-v2606"
BASE_CASE="$OPENFOAM_DIR/tutorials/incompressible/icoFoam/cavity/cavity"
PETSC_DIR="$THIRDPARTY_DIR/platforms/linux64IcxDPInt32/petsc-3.25.1"
PETSC_FOAM_LIB="$OPENFOAM_DIR/platforms/linux64IcxDPInt32Opt/lib/libpetscFoam.so"
NPROCS="${PBS_NCPUS:-8}"
MESH_N="${MESH_N:-256}"
REPEATS="${REPEATS:-3}"
DELTA_T="${DELTA_T:-0.00001}"
END_TIME="${END_TIME:-0.005}"
WRITE_INTERVAL="${WRITE_INTERVAL:-500}"
RUN_DIR="$OPENFOAM_ROOT/profiles/petsc_mkl_${PBS_JOBID:-manual}"
LOG_DIR="$OPENFOAM_ROOT/logs/profiles"

mkdir -p "$RUN_DIR" "$LOG_DIR"
exec > >(tee "$LOG_DIR/petsc_mkl_${PBS_JOBID:-manual}.out") 2>&1
trap 'rc=$?; echo "ERROR: line $LINENO: $BASH_COMMAND (exit $rc)" >&2' ERR

[[ "$NPROCS" == 8 ]] || { echo "ERROR: this configuration requires 8 MPI ranks, got $NPROCS"; exit 1; }
[[ -d "$BASE_CASE" ]] || { echo "ERROR: base case not found: $BASE_CASE"; exit 1; }
[[ -f "$PETSC_DIR/lib/libpetsc.so" ]] || { echo "ERROR: PETSc not found: $PETSC_DIR"; exit 1; }
[[ -f "$PETSC_FOAM_LIB" ]] || { echo "ERROR: petsc4Foam not found: $PETSC_FOAM_LIB"; exit 1; }

module purge
module load intel-compiler-llvm/2025.3.2
module load intel-mpi/2021.17.2
module load metis/5.1.0
module load intel-mkl/2025.3.1
module load boost/1.84.0
module load ipm/2.0.6

export WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx
export WM_MPLIB=INTELMPI I_MPI_CC=icx I_MPI_CXX=icpx
source "$OPENFOAM_DIR/etc/bashrc" || true

[[ "${WM_PROJECT_DIR:-}" == "$OPENFOAM_DIR" ]] || { echo "ERROR: OpenFOAM environment failed"; exit 1; }
[[ "$WM_OPTIONS" == linux64IcxDPInt32Opt ]] || { echo "ERROR: unexpected WM_OPTIONS: $WM_OPTIONS"; exit 1; }

export FOAM_MPI=intelmpi MPI_ARCH_PATH="$I_MPI_ROOT"
export PETSC_DIR PETSC_ARCH_PATH="$PETSC_DIR" PETSC_ARCH=""
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_THREADING_LAYER=SEQUENTIAL
export LIBRARY_PATH="$PETSC_DIR/lib:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="$FOAM_LIBBIN/$FOAM_MPI:$FOAM_LIBBIN:$PETSC_DIR/lib:$MKLROOT/lib:${LD_LIBRARY_PATH:-}"
export PKG_CONFIG_PATH="$PETSC_DIR/lib/pkgconfig:${PKG_CONFIG_PATH:-}"

IPM_LIBRARY=/apps/ipm/2.0.6/lib/mpich/libipm.so.0.0.0
[[ -f "$IPM_LIBRARY" ]] || { echo "ERROR: IPM library not found: $IPM_LIBRARY"; exit 1; }

PETSC_FOAM_LDD="$(ldd "$PETSC_FOAM_LIB")"
! grep -q 'not found' <<< "$PETSC_FOAM_LDD" || { echo "$PETSC_FOAM_LDD"; exit 1; }
for dependency in libpetsc libmkl_core libmkl_intel_lp64 libmkl_sequential libmpi; do
    grep -q "$dependency" <<< "$PETSC_FOAM_LDD" || {
        echo "ERROR: libpetscFoam dependency missing: $dependency"
        exit 1
    }
done

prepare_case()
{
    local case_dir="$1"
    local pressure_backend="$2"

    cp -a "$BASE_CASE" "$case_dir"
    sed -i "s/(20 20 1)/("$MESH_N" "$MESH_N" 1)/" "$case_dir/system/blockMeshDict"

    foamDictionary -entry numberOfSubdomains -set "$NPROCS" "$case_dir/system/decomposeParDict"
    foamDictionary -entry method -set hierarchical "$case_dir/system/decomposeParDict"
    foamDictionary -entry coeffs -remove "$case_dir/system/decomposeParDict"
    foamDictionary -entry hierarchicalCoeffs -set '{ n (4 2 1); order xyz; }' "$case_dir/system/decomposeParDict"
    foamDictionary -entry deltaT -set "$DELTA_T" "$case_dir/system/controlDict"
    foamDictionary -entry endTime -set "$END_TIME" "$case_dir/system/controlDict"
    foamDictionary -entry writeInterval -set "$WRITE_INTERVAL" "$case_dir/system/controlDict"

    foamDictionary -entry solvers/U/solver -set PBiCGStab "$case_dir/system/fvSolution"
    foamDictionary -entry solvers/U/smoother -remove "$case_dir/system/fvSolution"
    foamDictionary -entry solvers/U/preconditioner -set distributedDILU "$case_dir/system/fvSolution"

    if [[ "$pressure_backend" == native ]]; then
        foamDictionary -entry solvers/p/solver -set GAMG "$case_dir/system/fvSolution"
        foamDictionary -entry solvers/p/preconditioner -remove "$case_dir/system/fvSolution"
        foamDictionary -entry solvers/p/smoother -set GaussSeidel "$case_dir/system/fvSolution"
        foamDictionary -entry solvers/pFinal/solver -set GAMG "$case_dir/system/fvSolution"
        foamDictionary -entry solvers/pFinal/preconditioner -remove "$case_dir/system/fvSolution"
        foamDictionary -entry solvers/pFinal/smoother -set GaussSeidel "$case_dir/system/fvSolution"
    else
        foamDictionary -entry libs -set '(petscFoam)' "$case_dir/system/controlDict"
        foamDictionary -entry solvers/p -set \
            '{ solver petsc; tolerance 1e-06; relTol 0.05; petsc { verbose true; options { mat_type aijmkl; ksp_type cg; ksp_cg_single_reduction true; pc_type hypre; pc_hypre_type boomeramg; } caching { matrix { update always; } preconditioner { update always; } } } }' \
            "$case_dir/system/fvSolution"
        foamDictionary -entry solvers/pFinal -set \
            '{ solver petsc; tolerance 1e-06; relTol 0; petsc { verbose true; options { mat_type aijmkl; ksp_type cg; ksp_cg_single_reduction true; pc_type hypre; pc_hypre_type boomeramg; } caching { matrix { update always; } preconditioner { update always; } } } }' \
            "$case_dir/system/fvSolution"
    fi

    blockMesh -case "$case_dir" > "$case_dir/log.blockMesh" 2>&1
    decomposePar -case "$case_dir" > "$case_dir/log.decomposePar" 2>&1
}

check_log()
{
    local log="$1"
    grep -q '^End$' "$log" || { echo "ERROR: successful termination not found in $log"; exit 1; }
    ! grep -Eq 'FOAM FATAL|PETSC ERROR|Segmentation fault|KILLED BY SIGNAL: 8|Floating point exception \(core dumped\)' "$log" || {
        echo "ERROR: solver failure found in $log"
        exit 1
    }
}

run_repeats()
{
    local label="$1"
    local case_dir="$2"
    local require_petsc="$3"
    local mpi_exec ipm_dir ipm_xml

    mpi_exec="$(command -v mpiexec)"
    for ((run=1; run<=REPEATS; run++)); do
        echo "$label: timed run $run/$REPEATS"
        ipm_dir="$case_dir/ipm.$run"
        mkdir -p "$ipm_dir"
        /usr/bin/time \
            -f 'elapsed_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kb=%M' \
            -o "$case_dir/time.$run.txt" \
            "$mpi_exec" -n "$NPROCS" env \
                IPM_LOG=full IPM_REPORT=full IPM_LOGDIR="$ipm_dir" \
                LD_PRELOAD="$IPM_LIBRARY" \
                icoFoam -parallel -case "$case_dir" \
            > "$case_dir/log.icoFoam.$run" 2>&1
        check_log "$case_dir/log.icoFoam.$run"
        if [[ "$require_petsc" == true ]]; then
            grep -q 'PETSc-' "$case_dir/log.icoFoam.$run" || {
                echo "ERROR: PETSc solver was not observed in $case_dir/log.icoFoam.$run"
                exit 1
            }
            grep -q -- '-eqn_p_mat_type aijmkl' "$case_dir/log.icoFoam.$run" || {
                echo "ERROR: PETSc AIJMKL matrix option was not observed"
                exit 1
            }
        fi
        ipm_xml="$(find "$ipm_dir" -type f -name '*.ipm.xml' -print -quit)"
        [[ -n "$ipm_xml" ]] || { echo "ERROR: no IPM XML for $label run $run"; exit 1; }
        ipm_parse -full "$ipm_xml" > "$case_dir/ipm.$run.txt"
    done
}

prepare_case "$RUN_DIR/native_gamg" native
prepare_case "$RUN_DIR/petsc_mkl" petsc

run_repeats native_gamg "$RUN_DIR/native_gamg" false
run_repeats petsc_mkl "$RUN_DIR/petsc_mkl" true

average_elapsed()
{
    awk -F= '/^elapsed_seconds=/{sum += $2; count++} END {if (count) printf "%.6f", sum/count}' "$1"/time.*.txt
}

NATIVE_AVG="$(average_elapsed "$RUN_DIR/native_gamg")"
PETSC_AVG="$(average_elapsed "$RUN_DIR/petsc_mkl")"
RATIO="$(awk -v native="$NATIVE_AVG" -v petsc="$PETSC_AVG" 'BEGIN {if (petsc > 0) printf "%.4f", native/petsc}')"

{
    printf 'backend\taverage_seconds\n'
    printf 'native_gamg\t%s\n' "$NATIVE_AVG"
    printf 'petsc_mkl_aij\t%s\n' "$PETSC_AVG"
    printf 'native_over_petsc\t%s\n' "$RATIO"
} > "$RUN_DIR/summary.tsv"

echo
echo "================ PETSc/MKL benchmark summary ================"
column -t -s $'\t' "$RUN_DIR/summary.tsv" 2>/dev/null || cat "$RUN_DIR/summary.tsv"
echo "Case: icoFoam cavity, ${MESH_N}x${MESH_N}x1, $NPROCS MPI ranks, 500 steps"
echo "PETSc matrix type: AIJMKL"
echo "Results: $RUN_DIR"
