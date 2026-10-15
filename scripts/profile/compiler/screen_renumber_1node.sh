#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=00:30:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N renum_screen_1n
#PBS -j oe

# Screen renumberMesh methods on the prepared 104-rank occDrivAerStaticMesh case with
# renumberMesh -dry-run -frontWidth (computes the ordering, writes nothing). The case is already
# renumbered with the default CuthillMcKee, so "Before" = current ordering and each method's
# "After" = what it would give. Lower band / profile / rms frontwidth = neighbouring cells closer
# in memory. Methods: CuthillMcKee (current), RCM, spring, Sloan (libSloanRenumber).

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
CASE=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/renumber-screen/r1n_${PBS_JOBID%%.*}
fh=(-fileHandler hostCollated)

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "SCREEN_FAILED line=$LINENO status=$?"' ERR
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/../logFiles/$l.log"; done

for m in CuthillMcKee RCM spring Sloan; do
    {
        echo 'FoamFile { version 2.0; format ascii; class dictionary; object renumberMeshDict; }'
        if [[ $m == Sloan ]]; then echo 'libs (SloanRenumber);'; fi
        echo "method $m;"
        # spring requires its coefficients (annotated renumberMeshDict defaults)
        [[ $m == spring ]] && echo 'springCoeffs { maxCo 0.01; freezeFraction 0.999; maxIter 1000; }'
        true
    } > "$ROOT/renumberMeshDict.$m"
done

set +eu
module purge
module load OpenFOAM/2512
set -eu
cd "$CASE"
for m in CuthillMcKee RCM spring Sloan; do
    echo "DRYRUN_START $m utc=$(date -u +%FT%TZ)"
    status=0
    timeout 600 mpirun -np 104 --map-by core --bind-to core renumberMesh -constant -parallel "${fh[@]}" \
        -dry-run -frontWidth -no-fields -dict "$ROOT/renumberMeshDict.$m" > "$ROOT/dryrun_$m.log" 2>&1 || status=$?
    if [[ $status != 0 ]] || ! grep -q '^End' "$ROOT/dryrun_$m.log"; then
        echo "DRYRUN_FAILED $m status=$status: $(grep -m1 -iE 'FATAL|error|Unknown' "$ROOT/dryrun_$m.log" | cut -c1-160)"
        continue
    fi
    echo "DRYRUN $m method=$(grep -m1 -oE 'renumber-method: [A-Za-z]+' "$ROOT/dryrun_$m.log") $(awk '/^After renumbering/{f=1} f && /band|profile|frontwidth/ {printf "%s=%s ", $1 ($2=="frontwidth"?"_fw":""), $NF}' "$ROOT/dryrun_$m.log")"
done
cd "$ROOT"
[[ $(find "$CASE" -maxdepth 1 -newer "$ROOT/renumberMeshDict.Sloan" -name 'processors*' | wc -l) == 0 ]]   # dry-run wrote nothing

python3 - "$ROOT" <<'EOF'
import re, sys
from pathlib import Path
root = Path(sys.argv[1])
def stats(text, which):
    blk = text.split(which, 1)[1] if which in text else ''
    get = lambda k: (re.search(rf'{k}\s*:\s*([0-9.eE+-]+)', blk) or [None, None])[1]
    return get('band'), get('profile'), get('rms frontwidth')
lines = ['# renumberMesh method screen, occDrivAerStaticMesh 65M, 104-rank scotch case (dry-run)', '',
         '| Method | band | profile | rms frontwidth |', '|---|---:|---:|---:|']
before = None
for m in ('CuthillMcKee', 'RCM', 'spring', 'Sloan'):
    f = root / f'dryrun_{m}.log'
    t = f.read_text() if f.exists() else ''
    if before is None and 'Before renumbering' in t: before = stats(t, 'Before renumbering')
    b, p, w = stats(t, 'After renumbering')
    lines.append(f'| {m} | {b} | {p} | {w} |' if b else f'| {m} | failed | | |')
if before: lines.insert(4, f'| current ordering (CuthillMcKee from prep) | {before[0]} | {before[1]} | {before[2]} |')
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "SCREEN_COMPLETE root=$ROOT"
