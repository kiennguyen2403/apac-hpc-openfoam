#!/bin/bash -l
#PBS -N openmpi6_build
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=32
#PBS -l mem=64GB
#PBS -l walltime=01:00:00
#PBS -l jobfs=12GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -j oe
#PBS -o /scratch/ya41/kn3896/openfoam/logs/build

set -eo pipefail
: "${PBS_JOBFS:?Submit this script with qsub}"

module purge
module load intel-compiler-llvm/2025.3.2 autoconf/2.72 python3/3.11.7
module load ucx/1.18.1 hwloc/2.12.1 libevent/2.1.12 pmix/6.0.0 prrte/4.0.0

src=/home/565/kn3896/src/openmpi-6.0.0rc2
prefix=/scratch/ya41/kn3896/openfoam/ThirdParty-v2512/platforms/linux64Icx/openmpi-6.0.0rc2
build_dir="$PBS_JOBFS/openmpi6"
log_dir="/scratch/ya41/kn3896/openfoam/logs/build/openmpi6_${PBS_JOBID}"

if [[ -e "$prefix" ]]; then
    echo "Install directory already exists: $prefix" >&2
    exit 1
fi
mkdir -p "$log_dir"
echo "Build logs: $log_dir"

export CC=icx CXX=icpx
export CFLAGS="-Ofast -march=sapphirerapids -mprefer-vector-width=512 -fPIC"
export CXXFLAGS="$CFLAGS"
cp -a "$src" "$build_dir"
cd "$build_dir"

./autogen.pl --no-3rdparty=libevent,hwloc,pmix,prrte > "$log_dir/autogen.log" 2>&1
./configure --prefix="$prefix" \
    --with-ucx=/apps/ucx/1.18.1 \
    --with-hwloc=/apps/hwloc/2.12.1 \
    --with-libevent=/apps/libevent/2.1.12 \
    --with-pmix=/apps/pmix/6.0.0 \
    --with-prrte=/apps/prrte/4.0.0 \
    --disable-mpi-fortran --disable-oshmem \
    --disable-debug --disable-static > "$log_dir/configure.log" 2>&1
make -j "$PBS_NCPUS" > "$log_dir/make.log" 2>&1
make install > "$log_dir/install.log" 2>&1

"$prefix/bin/ompi_info" --version
"$prefix/bin/ompi_info" --param coll tuned --level 9 > "$log_dir/collectives.log"
grep -F 'allgather_reduce' "$log_dir/collectives.log"
"$prefix/bin/ompi_info" --param pml ucx > "$log_dir/ucx.log"
grep -F 'MCA pml: ucx' "$log_dir/ucx.log"
prte_info > "$log_dir/prrte.log"
grep -F 'plm: tm' "$log_dir/prrte.log"
echo "Installed Open MPI 6 in $prefix"
