#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N renum_raw_screen
#PBS -j oe

# Option A, step 2: renumberMesh -dry-run -frontWidth on the RAW 104-rank decomposition
# (case104_raw), i.e. a single pass as in the mandated chain. Question: does spring alone reach
# the ~2.3e11 profile it gave on top of CuthillMcKee (job 180787768), or does it need CM first?
# Methods: CuthillMcKee (team default), spring. Writes nothing.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
C=/scratch/ya41/kn3896/openfoam/contest-results/renumber-screen/case104_raw
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/renumber-screen/raw_${PBS_JOBID%%.*}
mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "SCREEN_FAILED line=$LINENO status=$?"' ERR
grep -q '^End' "$C/logFiles/10_decomposePar.log"

printf 'FoamFile { version 2.0; format ascii; class dictionary; object renumberMeshDict; }\nmethod CuthillMcKee;\n' > "$ROOT/renumberMeshDict.CuthillMcKee"
printf 'FoamFile { version 2.0; format ascii; class dictionary; object renumberMeshDict; }\nmethod spring;\nspringCoeffs { maxCo 0.01; freezeFraction 0.999; maxIter 1000; }\n' > "$ROOT/renumberMeshDict.spring"

set +eu
module purge
module load OpenFOAM/2512
set -eu
cd "$C"
for m in CuthillMcKee spring; do
    echo "DRYRUN_START $m utc=$(date -u +%FT%TZ)"
    status=0
    timeout 1200 mpirun -np 104 --map-by core --bind-to core renumberMesh -constant -parallel -fileHandler collated \
        -dry-run -frontWidth -no-fields -dict "$ROOT/renumberMeshDict.$m" > "$ROOT/dryrun_$m.log" 2>&1 || status=$?
    if [[ $status != 0 ]] || ! grep -q '^End' "$ROOT/dryrun_$m.log"; then
        echo "DRYRUN_FAILED $m status=$status: $(grep -m1 -iE 'FATAL|error' "$ROOT/dryrun_$m.log" | cut -c1-160)"; continue
    fi
    echo "DRYRUN $m $(awk '/^Before renumbering/{b=1;a=0} /^After renumbering/{a=1;b=0} (a||b) && /band|profile|frontwidth/ {printf "%s_%s=%s ", (a?"after":"before"), $1, $NF}' "$ROOT/dryrun_$m.log")"
done
cd "$ROOT"
{
    echo '# renumberMesh single pass on the RAW 104-rank decomposition (dry-run)'
    echo
    echo 'Reference: on top of CuthillMcKee (job 180787768) spring gave profile 2.30e11, rms frontwidth 3890; CuthillMcKee 5.44e11 / 9491.'
    echo
    for m in CuthillMcKee spring; do echo "$m: $(grep "^DRYRUN $m" job.log | cut -d' ' -f3-)"; done
} > REPORT.md
cat REPORT.md
echo "SCREEN_COMPLETE root=$ROOT"
