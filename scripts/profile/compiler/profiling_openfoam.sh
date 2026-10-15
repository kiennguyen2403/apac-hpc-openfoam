#!/bin/bash -l
#PBS -N compiler_smoke
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104
#PBS -l mem=64GB
#PBS -l walltime=00:10:00
#PBS -l jobfs=4GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -j oe

# One-node ICX / default OpenFOAM module / Clang DrivAer comparison.
# Default: 6-step tiny smoke. Full: -v BENCH_MODE=full,SMOKE_EVIDENCE,SMOKE_PBS_OUTPUT
# and -l mem=400GB,walltime=01:00:00; require afterok dependency on smoke.
# Both modes use all 104 physical cores, two repeats and per-core CPU samples.
set -euo pipefail
export BENCH_MODE=${BENCH_MODE:-smoke}
exec bash -l /home/565/kn3896/drivaer-contest-profile/compiler-one-node/run-benchmark.sh
