#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=416,mem=400gb
#PBS -l walltime=00:20:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N sloan_dry_4n
#PBS -j oe

# BONUS screen: renumberMesh -dry-run with Sloan (third-party libSloanRenumber, shipped with the
# site module) vs the default CuthillMcKee on the prepared 416-part scotch (4 104)/(1 0.5) case.
# Sloan crashed on the 104-part case (rank 85: "Map is not one-to-one"); the 416 split has
# different subdomains. Nothing is written (-dry-run). The case is already CuthillMcKee-ordered
# (renumbered in place by the prep chain), so both methods start from that ordering.

set -uo pipefail
: "${PBS_JOBID:?submit through PBS}"
C=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/case416_w05
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/renumber-screen/sloan416_${PBS_JOBID%%.*}
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
[[ $(sort -u "$PBS_NODEFILE" | wc -l) == 4 ]] || exit 1
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$C/logFiles/$l.log" || exit 1; done
printf 'FoamFile { version 2.0; format ascii; class dictionary; object renumberMeshDict; }\nmethod CuthillMcKee;\n' > "$ROOT/renumberMeshDict.CuthillMcKee"
printf 'FoamFile { version 2.0; format ascii; class dictionary; object renumberMeshDict; }\nlibs (SloanRenumber);\nmethod Sloan;\n' > "$ROOT/renumberMeshDict.Sloan"

module purge
module load OpenFOAM/2512
cd "$C"
for m in CuthillMcKee Sloan; do
    echo "DRYRUN_START $m utc=$(date -u +%FT%TZ)"
    status=0
    timeout 600 mpirun -np 416 --map-by ppr:104:node --bind-to core renumberMesh -constant -parallel "${fh[@]}" \
        -dry-run -frontWidth -no-fields -dict "$ROOT/renumberMeshDict.$m" > "$ROOT/dryrun_$m.log" 2>&1 || status=$?
    if [[ $status != 0 ]] || ! grep -q '^End' "$ROOT/dryrun_$m.log"; then
        echo "DRYRUN_FAILED $m status=$status: $(grep -m1 -E 'FOAM FATAL|Map is not|error' "$ROOT/dryrun_$m.log" | cut -c1-160)"
        grep -m3 -B1 -A3 'FOAM FATAL' "$ROOT/dryrun_$m.log" | head -12
        continue
    fi
    echo "DRYRUN $m $(grep -m1 -oE 'renumber-method: [A-Za-z]+' "$ROOT/dryrun_$m.log") $(awk '/^Before renumbering/{b=1;a=0} /^After renumbering/{a=1;b=0} (a||b) && /band|profile|frontwidth/ {printf "%s_%s=%s ", (a?"after":"before"), $1, $NF}' "$ROOT/dryrun_$m.log")"
done
echo "SLOAN_DRYRUN_COMPLETE root=$ROOT"
