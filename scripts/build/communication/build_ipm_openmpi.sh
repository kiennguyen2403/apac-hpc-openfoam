#!/bin/bash -l
#PBS -N ipm_openmpi_fix
#PBS -P ya41
#PBS -q normal
#PBS -l ncpus=2
#PBS -l mem=8GB
#PBS -l walltime=00:10:00
#PBS -l jobfs=2GB
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -j oe

set -eEo pipefail
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR
: "${PBS_JOBID:?Submit this build and validation with qsub}"
: "${PBS_JOBFS:?PBS jobfs required}"
[[ $PBS_NCPUS == 2 ]]

IPM_SOURCE=/scratch/ya41/kn3896/openfoam/ThirdParty-v2512/ipm-2.0.6-openmpi
IPM_PLATFORM=/scratch/ya41/kn3896/openfoam/ThirdParty-v2512/platforms/linux64Gcc
IPM_BUILD_LOG=/scratch/ya41/kn3896/openfoam/logs/build/ipm_openmpi_$PBS_JOBID
PROFILE_ENTRY=/scratch/ya41/kn3896/openfoam/scripts/profile/communication/profiling_openfoam.sh
PROFILE_RESULTS=/scratch/ya41/kn3896/openfoam/contest-results/communication
mkdir -p "$IPM_BUILD_LOG"
[[ $(git -C "$IPM_SOURCE" rev-parse HEAD) == b008141ee16d39b33e20bffde615564afa107575 ]]

module purge
module load gcc/11.1.0 openmpi/5.0.8
unset LD_PRELOAD
# NCI compiler wrappers supply language-specific MPI/GNU include and library paths.
export OMPI_CC=gcc OMPI_CXX=g++ OMPI_FC=gfortran
# IPM skips compiler detection when MPICC is preset, leaving libtool flags unset.
unset MPICC MPIF77 MPIFC MPILIBS
export CC=gcc CXX=g++ F77=gfortran FC=gfortran
export CFLAGS='-O2 -g -fcommon'
module list > "$IPM_BUILD_LOG/modules.txt" 2>&1
gcc --version > "$IPM_BUILD_LOG/compiler.txt"
git -C "$IPM_SOURCE" status --short > "$IPM_BUILD_LOG/source-status-before.txt"
(cd "$IPM_SOURCE"; autoreconf --install --force) > "$IPM_BUILD_LOG/bootstrap.log" 2>&1

# Both copies use the same source, compiler, MPI and PAPI. Only the upstream
# HAVE_NO_DELAYED_FINALIZE build switch differs, to isolate shutdown timing.
for finalize_mode in delayed eager; do
    install_dir=$IPM_PLATFORM/ipm-2.0.6-openmpi-5.0.8-mpi3-$finalize_mode
    build_dir=$PBS_JOBFS/ipm-build-$finalize_mode
    [[ ! -e $install_dir && ! -e $build_dir ]] || {
        echo "ERROR: refusing to overwrite an existing IPM build: $install_dir" >&2
        exit 2
    }
    mkdir -p "$build_dir"
    export BUILD_ROOT=$build_dir
    export CPPFLAGS=
    if [[ $finalize_mode == eager ]]; then export CPPFLAGS=-DHAVE_NO_DELAYED_FINALIZE=1; fi
    echo "IPM_BUILD_START mode=$finalize_mode"
    (cd "$build_dir"; "$IPM_SOURCE/configure" --prefix="$install_dir" \
        --enable-shared --disable-static --disable-parser --with-papi=/apps/papi/5.7.0) \
        > "$IPM_BUILD_LOG/configure-$finalize_mode.log" 2>&1
    # Preserve the site module\'s MPI-3 wrappers when generating the profiler.
    make -C "$build_dir" -j 2 IPM_KEYFILE_MPI=/apps/ipm/2.0.6/etc/ipm_key_mpi \
        > "$IPM_BUILD_LOG/make-$finalize_mode.log" 2>&1
    make -C "$build_dir" install IPM_KEYFILE_MPI=/apps/ipm/2.0.6/etc/ipm_key_mpi \
        > "$IPM_BUILD_LOG/install-$finalize_mode.log" 2>&1
    install -m 644 /apps/ipm/2.0.6/etc/ipm_key_mpi "$install_dir/etc/ipm_key_mpi"
    nm -D "$install_dir/lib/libipm.so" > "$IPM_BUILD_LOG/symbols-$finalize_mode.txt"
    grep -Eq " T MPI_Iallreduce$" "$IPM_BUILD_LOG/symbols-$finalize_mode.txt"
    grep -Eq " T MPI_Ibcast$" "$IPM_BUILD_LOG/symbols-$finalize_mode.txt"
    [[ -s $install_dir/lib/libipm.so ]]
    {
        echo 'source_commit=b008141ee16d39b33e20bffde615564afa107575'
        echo "finalize_mode=$finalize_mode"
        echo "cppflags=$CPPFLAGS"
        echo "build_job=$PBS_JOBID"
        sha256sum "$install_dir/lib/libipm.so"
        readelf -d "$install_dir/lib/libipm.so"
    } > "$IPM_BUILD_LOG/library-$finalize_mode.txt"
    echo "IPM_BUILD_SUCCESS mode=$finalize_mode prefix=$install_dir"
done

# First test the otherwise-identical delayed build. Preserve failures as
# diagnostic evidence; only the eager build must pass to complete this job.
if env PROFILE_MODE=ipm-smoke PROFILE_VARIANTS=openmpi \
    RESULT_ROOT="$PROFILE_RESULTS/ipm_openmpi_delayed_control" \
    IPM_OPENMPI_LIBRARY="$IPM_PLATFORM/ipm-2.0.6-openmpi-5.0.8-mpi3-delayed/lib/libipm.so" \
    "$PROFILE_ENTRY" > "$IPM_BUILD_LOG/delayed-control.log" 2>&1; then
    delayed_status=0
else
    delayed_status=$?
fi
echo "DELAYED_CONTROL_STATUS=$delayed_status"

# Exercise retained MPI-3 wrappers and validate their exact recorded counts.
collective_dir=$PROFILE_RESULTS/ipm_collectives_$PBS_JOBID
mkdir -p "$collective_dir"
mpicc /home/565/kn3896/drivaer-contest-profile/tests/ipm-collectives.c -o "$PBS_JOBFS/ipm-collectives"
module load ipm/2.0.6
unset LD_PRELOAD
export IPM_LOG=full IPM_REPORT=full IPM_LOGWRITER=serial IPM_LOGDIR=$collective_dir
export IPM_KEYFILE=$IPM_PLATFORM/ipm-2.0.6-openmpi-5.0.8-mpi3-eager/etc/ipm_key_mpi
timeout --kill-after=5s 30s mpirun -np 2 --map-by core --bind-to core \
    -x LD_PRELOAD="$IPM_PLATFORM/ipm-2.0.6-openmpi-5.0.8-mpi3-eager/lib/libipm.so" \
    "$PBS_JOBFS/ipm-collectives" > "$collective_dir/test.log" 2>&1
collective_xml=("$collective_dir"/*.xml)
[[ ${#collective_xml[@]} == 1 && -s ${collective_xml[0]} ]]
ipm_parse -full "${collective_xml[0]}" > "$collective_dir/report.txt" 2>&1
grep -Eq "MPI_Iallreduce[[:space:]]+[^[:space:]]+[[:space:]]+16[[:space:]]" "$collective_dir/report.txt"
grep -Eq "MPI_Ibcast[[:space:]]+[^[:space:]]+[[:space:]]+16[[:space:]]" "$collective_dir/report.txt"
echo "MPI3_WRAPPER_TEST_SUCCESS results=$collective_dir"

# The normal helper selects the eager copy by default. Give the second test
# independent working and output directories without deleting control inputs.
mkdir "$PBS_JOBFS/fixed-mock"
env PBS_JOBFS="$PBS_JOBFS/fixed-mock" PROFILE_MODE=ipm-smoke \
    PROFILE_VARIANTS='intelmpi openmpi' RESULT_ROOT="$PROFILE_RESULTS" \
    "$PROFILE_ENTRY" > "$IPM_BUILD_LOG/fixed-mock.log" 2>&1
echo "IPM_OPENMPI_FIX_SUCCESS results=$PROFILE_RESULTS/ipm_smoke_$PBS_JOBID logs=$IPM_BUILD_LOG"
