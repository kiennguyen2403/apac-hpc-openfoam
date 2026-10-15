#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=26,mem=128gb
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N dw_screen_416
#PBS -j oe

# Screen scotch domainWeights for the 4-node (416 = 4 x 104) decomposition of
# occDrivAerStaticMesh 65M with decomposePar -dry-run (computes the partition, writes nothing).
# OpenFOAM turns domainWeights into integer SCOTCH_archTleaf link costs (w/min(w)): the team's
# (1 0.01) tells scotch cross-node edges cost 100x within-node ones, so within-node
# communication balance is barely optimised (worst rank 58 neighbours vs median ~11).
# Variants: (1 0.01) current, (1 0.2), (1 0.5), (1 1), and plain scotch (no domains).
# Reported per variant: cells / neighbours / processor faces min, max, median (max/median %).

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/dw416_${PBS_JOBID%%.*}
C=$ROOT/case

mkdir -p "$C/constant"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "SCREEN_FAILED line=$LINENO status=$?"' ERR
cp -a "$CASE_SRC/system" "$C/"
cp -a "$CASE_SRC"/constant/{transportProperties,turbulenceProperties} "$C/constant/"
ln -s "$MESH" "$C/constant/polyMesh"
cp "$TEAM/dicts/4node/caseDefinition"   "$C/system/include/caseDefinition"
cp "$TEAM/dicts/4node/decomposeParDict" "$C/system/decomposeParDict"
grep -E '^(nCores|decompositionMethod|nDomains|nDomainWeights)' "$C/system/include/caseDefinition"

set +eu
module purge
module load OpenFOAM/2512
set -eu
cd "$C"
cp system/include/caseDefinition system/include/caseDefinition.team
cp system/decomposeParDict system/decomposeParDict.team

for v in w0.01 w0.2 w0.5 w1 plain; do
    cp system/include/caseDefinition.team system/include/caseDefinition
    cp system/decomposeParDict.team system/decomposeParDict
    case $v in
        w*)    sed -i "s/^nDomainWeights .*/nDomainWeights (1 ${v#w});/" system/include/caseDefinition ;;
        # plain scotch: drop the multi-level coefficients
        plain) sed -i '/^scotchCoeffs/,/^}/{/domains\|domainWeights/d}' system/decomposeParDict ;;
    esac
    {
        echo "== variant $v"
        foamDictionary -expand system/decomposeParDict | sed -n '/^numberOfSubdomains/p;/^method/p;/^scotchCoeffs/,/^}/p'
    } > "$ROOT/dict_$v.txt"
    [[ $(awk '/^numberOfSubdomains/ {sub(";","",$2); print $2}' "$ROOT/dict_$v.txt") == 416 ]]
    echo "DRYRUN_START $v utc=$(date -u +%FT%TZ)"
    decomposePar -dry-run > "$ROOT/dryrun_$v.log" 2>&1
    grep -q '^End' "$ROOT/dryrun_$v.log"
    echo "DRYRUN $v $(grep -E '^(Cells|Neigh|Faces) ' "$ROOT/dryrun_$v.log" | tr '\n' ' ')"
done
cp system/include/caseDefinition.team system/include/caseDefinition
cp system/decomposeParDict.team system/decomposeParDict

python3 - "$ROOT" <<'EOF'
import re, sys
from pathlib import Path
root = Path(sys.argv[1])
pat = re.compile(r'^(Cells|Neigh|Faces) min:(\d+) max:(\d+) median:(\d+)', re.M)
lines = ['# scotch domainWeights screen, 416 = 4 x 104, occDrivAerStaticMesh 65M (decomposePar -dry-run)', '',
         '| Variant | Cells max/median | Neigh min / median / max | Faces min / median / max | Faces max/median |',
         '|---|---:|---:|---:|---:|']
for v in ('w0.01', 'w0.2', 'w0.5', 'w1', 'plain'):
    s = {m[0]: tuple(map(int, m[1:])) for m in pat.findall((root / f'dryrun_{v}.log').read_text())}
    c, n, f = s['Cells'], s['Neigh'], s['Faces']
    lines.append(f'| {v} | {100*c[1]/c[2]:.2f}% | {n[0]} / {n[2]} / **{n[1]}** | {f[0]} / {f[2]} / **{f[1]}** | {100*f[1]/f[2]:.1f}% |')
lines += ['', 'w0.01 = team setting (domains (4 104), domainWeights (1 0.01)). Lower worst-case Neigh/Faces = less MPI imbalance.',
          'Not shown: how many faces cross nodes (a lower weight trades more cross-node faces for better within-node balance).']
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "SCREEN_COMPLETE root=$ROOT"
