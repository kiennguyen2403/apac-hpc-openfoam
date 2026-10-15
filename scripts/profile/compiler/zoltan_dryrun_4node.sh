#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:20:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N zoltan_dry_4n
#PBS -j oe

# BONUS screen: renumberMesh -dry-run with Zoltan (third-party; OpenFOAM zoltanRenumber, default
# ORDER_METHOD=LOCAL_HSFC = Hilbert space-filling curve through cell centres) vs the default
# CuthillMcKee, on the prepared 416-part scotch (4 104)/(1 0.5) case. Zoltan 3.901 built with
# ICX + Intel MPI (tools/zoltan-3.901-icx-impi); libzoltanRenumber in the ICX OptSPR tree, so both
# methods run with that renumberMesh + Intel MPI. Nothing is written. The case is already
# CuthillMcKee-ordered (renumbered in place by the prep chain).

set -uo pipefail
: "${PBS_JOBID:?submit through PBS}"
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
C=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/renumber-screen/zoltan416_${PBS_JOBID%%.*}
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]] || exit 1
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$C/logFiles/$l.log" || exit 1; done
printf 'FoamFile { version 2.0; format ascii; class dictionary; object renumberMeshDict; }\nmethod CuthillMcKee;\n' > "$ROOT/renumberMeshDict.CuthillMcKee"
printf 'FoamFile { version 2.0; format ascii; class dictionary; object renumberMeshDict; }\nlibs (zoltanRenumber);\nmethod zoltan;\nzoltanCoeffs { ORDER_METHOD LOCAL_HSFC; }\n' > "$ROOT/renumberMeshDict.zoltan"

set +u   # OpenFOAM bashrc is not nounset-safe
module purge
module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
export CC=icx CXX=icpx I_MPI_CC=icx I_MPI_CXX=icpx
source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=Icx WM_CC=icx WM_CXX=icpx \
    WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 PTSCOTCH_ARCH_PATH=none
set -u
[[ $WM_OPTIONS == linux64IcxDPInt32OptSPR ]] || exit 1
ls "$FOAM_LIBBIN/intelmpi/libzoltanRenumber.so" || exit 1
export I_MPI_LIBRARY_KIND=release I_MPI_FABRICS=shm:ofi I_MPI_OFI_PROVIDER=mlx I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0
export OMP_NUM_THREADS=1

cd "$C"
for m in CuthillMcKee zoltan; do
    echo "DRYRUN_START $m utc=$(date -u +%FT%TZ)"
    status=0
    timeout 600 mpirun -np 416 -ppn 104 renumberMesh -constant -parallel "${fh[@]}" \
        -dry-run -frontWidth -no-fields -dict "$ROOT/renumberMeshDict.$m" > "$ROOT/dryrun_$m.log" 2>&1 || status=$?
    if [[ $status != 0 ]] || ! grep -q '^End' "$ROOT/dryrun_$m.log"; then
        echo "DRYRUN_FAILED $m status=$status: $(grep -m1 -E 'FOAM FATAL|Zoltan|error' "$ROOT/dryrun_$m.log" | cut -c1-160)"
        grep -m2 -A4 'FOAM FATAL' "$ROOT/dryrun_$m.log" | head -12
        continue
    fi
    echo "DRYRUN $m $(grep -m1 -oE 'renumber-method: [A-Za-z]+' "$ROOT/dryrun_$m.log") $(grep -m1 -oE 'ORDER_METHOD = [A-Z_]+' "$ROOT/dryrun_$m.log") $(awk '/^Before renumbering/{b=1;a=0} /^After renumbering/{a=1;b=0} (a||b) && /band|profile|frontwidth/ {printf "%s_%s=%s ", (a?"after":"before"), $1, $NF}' "$ROOT/dryrun_$m.log")"
done
echo "ZOLTAN_DRYRUN_COMPLETE root=$ROOT"
