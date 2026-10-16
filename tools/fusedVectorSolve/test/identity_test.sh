#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=4,mem=16gb,walltime=00:15:00,storage=scratch/ya41
#PBS -N fused_identity
#PBS -j oe

# Bit-identity test of libfusedVectorSolve against stock OpenFOAM on the 3D box case
# (make_box_case.sh), serial and 4 ranks. usage: qsub -v WORK=<empty scratch dir>[,PRELOAD=<lib[:lib]>] identity_test.sh
# PRELOAD defaults to libfusedVectorSolve; any preloaded library must leave results bit-identical.
set -euo pipefail
: "${WORK:?qsub -v WORK=<dir>}"
T=/scratch/ya41/kn3896/openfoam/tools/fusedVectorSolve
FL=${PRELOAD:-$T/lib/linux64IcxDPInt32OptSPR/libfusedVectorSolve.so}
set +eu
module purge; module load intel-compiler-llvm/2025.3.2 intel-mpi/2021.17.2 > /dev/null 2>&1
source /scratch/ya41/kn3896/openfoam/OpenFOAM-v2512/etc/bashrc WM_COMPILER_TYPE=system WM_COMPILER=Icx \
    WM_CC=icx WM_CXX=icpx WM_MPLIB=INTELMPI WM_COMPILE_OPTION=OptSPR > /dev/null 2>&1
set -eu
export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512 I_MPI_FABRICS=shm FOAM_SIGFPE=true FUSED_VECTOR_SOLVE_REPORT=1 CSR_AMUL_REPORT=1
mkdir -p "$WORK" && cd "$WORK"
for l in ${FL//:/ }; do echo "library $l sha256 $(sha256sum "$l" | cut -d' ' -f1)"; done
for v in stock fused; do
    bash "$T/test/make_box_case.sh" "$WORK/$v" > /dev/null
    pre=""; [[ $v == fused ]] && pre=$FL
    (cd "$v" && blockMesh > log.blockMesh 2>&1 && cp -r . "../${v}_par"
     LD_PRELOAD=$pre simpleFoam > log.serial 2>&1
     echo "$v serial exit=$? steps=$(grep -c '^Time =' log.serial) $(grep -E 'fusedVectorSolve|csrAmul' log.serial | sort -u | tr '\n' ' ' || true)")
    (cd "${v}_par" && decomposePar -force > log.decompose 2>&1 && mpirun -np 4 -genv LD_PRELOAD "$pre" simpleFoam -parallel > log.par 2>&1
     echo "$v parallel exit=$? steps=$(grep -c '^Time =' log.par) $(grep -E 'fusedVectorSolve|csrAmul' log.par | sort -u | tr '\n' ' ' || true)"
     reconstructPar -latestTime > log.recon 2>&1)
done
for m in serial par; do
    s=stock; f=fused; [[ $m == par ]] && s=stock_par f=fused_par
    r=""; for fld in U p k omega nut; do cmp -s "$s/20/$fld" "$f/20/$fld" || r="$r $fld"; done
    if diff -q <(grep -E 'Solving for|continuity' "$s/log.$m") <(grep -E 'Solving for|continuity' "$f/log.$m") > /dev/null; then l=identical; else l=DIFFER; fi
    echo "IDENTITY $m: fields ${r:-all identical}; solver lines $l"
done
