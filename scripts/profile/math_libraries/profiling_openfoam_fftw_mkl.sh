#!/bin/bash -l
#PBS -N benchmark_fftw_mkl
#PBS -q normal
#PBS -l ncpus=1
#PBS -l mem=4GB
#PBS -l walltime=00:10:00
#PBS -l wd
#PBS -j oe
#PBS -o /home/565/kn3896/scratch/openfoam/logs/profiles/benchmark_fftw_mkl.out

set -eEo pipefail

OPENFOAM_ROOT=/home/565/kn3896/scratch/openfoam
SOURCE="$OPENFOAM_ROOT/benchmarks/fftw/fftw_execution_benchmark.c"
STANDARD_FFTW=/apps/fftw3/3.3.10
STANDARD_LIB="$STANDARD_FFTW/lib/Intel"
STANDARD_RUNTIME_LIB="$STANDARD_FFTW/lib"
MKL_FFTW="$OPENFOAM_ROOT/ThirdParty-v2606/platforms/linux64IcxDPInt32/fftw-mkl-2025.3.1"
GRID_SIZE="${GRID_SIZE:-256}"
REPETITIONS="${REPETITIONS:-20}"
RUN_DIR="$OPENFOAM_ROOT/profiles/fftw_mkl_${PBS_JOBID:-manual}"
LOG_DIR="$OPENFOAM_ROOT/logs/profiles"

mkdir -p "$RUN_DIR" "$LOG_DIR"
exec > >(tee "$LOG_DIR/fftw_mkl_${PBS_JOBID:-manual}.out") 2>&1
trap 'rc=$?; echo "ERROR: line $LINENO: $BASH_COMMAND (exit $rc)" >&2' ERR

module purge
module load intel-compiler-llvm/2025.3.2
module load intel-mkl/2025.3.1

[[ -f "$SOURCE" ]] || { echo "ERROR: source not found: $SOURCE"; exit 1; }
[[ -f "$STANDARD_LIB/libfftw3.so" ]] || { echo "ERROR: standard FFTW not found"; exit 1; }
[[ -f "$MKL_FFTW/lib/libfftw3.so" ]] || { echo "ERROR: MKL FFTW interface not found"; exit 1; }

icx -O3 -std=c11 -I"$STANDARD_FFTW/include" "$SOURCE" \
    -L"$STANDARD_LIB" -Wl,-rpath,"$STANDARD_RUNTIME_LIB" -lfftw3 -lm \
    -o "$RUN_DIR/benchmark_standard_fftw"

icx -O3 -std=c11 -I"$MKL_FFTW/include" "$SOURCE" \
    -L"$MKL_FFTW/lib" -Wl,-rpath,"$MKL_FFTW/lib" -lfftw3 -lm \
    -o "$RUN_DIR/benchmark_mkl_fftw"

STANDARD_LDD="$(LD_LIBRARY_PATH="$STANDARD_RUNTIME_LIB:${LD_LIBRARY_PATH:-}" ldd "$RUN_DIR/benchmark_standard_fftw")"
MKL_LDD="$(LD_LIBRARY_PATH="$MKL_FFTW/lib:$MKLROOT/lib:${LD_LIBRARY_PATH:-}" ldd "$RUN_DIR/benchmark_mkl_fftw")"

grep -q "$STANDARD_RUNTIME_LIB/libfftw3_Intel" <<< "$STANDARD_LDD" || { echo "$STANDARD_LDD"; exit 1; }
! grep -qi 'mkl' <<< "$STANDARD_LDD" || { echo "ERROR: standard FFTW unexpectedly resolves MKL"; exit 1; }
grep -q "$MKL_FFTW/lib/libfftw3" <<< "$MKL_LDD" || { echo "$MKL_LDD"; exit 1; }
grep -q 'libmkl_rt' <<< "$MKL_LDD" || { echo "ERROR: MKL FFTW interface does not resolve libmkl_rt"; exit 1; }
! grep -q 'not found' <<< "$STANDARD_LDD$MKL_LDD" || { echo "$STANDARD_LDD$MKL_LDD"; exit 1; }

export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_DYNAMIC=FALSE

echo "Running standard FFTW 3.3.10"
LD_LIBRARY_PATH="$STANDARD_RUNTIME_LIB:${LD_LIBRARY_PATH:-}" \
    "$RUN_DIR/benchmark_standard_fftw" "$GRID_SIZE" "$REPETITIONS" \
    | tee "$RUN_DIR/standard_fftw.txt"

echo "Running FFTW3 interface backed by Intel MKL"
LD_LIBRARY_PATH="$MKL_FFTW/lib:$MKLROOT/lib:${LD_LIBRARY_PATH:-}" \
    "$RUN_DIR/benchmark_mkl_fftw" "$GRID_SIZE" "$REPETITIONS" \
    | tee "$RUN_DIR/mkl_fftw.txt"

STANDARD_TIME="$(awk -F= '$1=="average_seconds"{print $2}' "$RUN_DIR/standard_fftw.txt")"
MKL_TIME="$(awk -F= '$1=="average_seconds"{print $2}' "$RUN_DIR/mkl_fftw.txt")"
STANDARD_CHECKSUM="$(awk -F= '$1=="checksum"{print $2}' "$RUN_DIR/standard_fftw.txt")"
MKL_CHECKSUM="$(awk -F= '$1=="checksum"{print $2}' "$RUN_DIR/mkl_fftw.txt")"
SPEEDUP="$(awk -v standard="$STANDARD_TIME" -v mkl="$MKL_TIME" 'BEGIN {printf "%.4f", standard/mkl}')"
RELATIVE_ERROR="$(awk -v standard="$STANDARD_CHECKSUM" -v mkl="$MKL_CHECKSUM" 'BEGIN {d=standard-mkl; if (d<0)d=-d; printf "%.3e", d/(standard?standard:1)}')"

{
    printf 'backend\taverage_seconds\n'
    printf 'standard_fftw\t%s\n' "$STANDARD_TIME"
    printf 'mkl_fftw_interface\t%s\n' "$MKL_TIME"
    printf 'mkl_speedup\t%s\n' "$SPEEDUP"
    printf 'checksum_relative_error\t%s\n' "$RELATIVE_ERROR"
} > "$RUN_DIR/summary.tsv"

echo
echo "================ FFTW/MKL benchmark summary ================"
column -t -s $'\t' "$RUN_DIR/summary.tsv" 2>/dev/null || cat "$RUN_DIR/summary.tsv"
echo "Transform: ${GRID_SIZE}x${GRID_SIZE}x${GRID_SIZE} real-to-complex, $REPETITIONS executions, one thread"
echo "Results: $RUN_DIR"
