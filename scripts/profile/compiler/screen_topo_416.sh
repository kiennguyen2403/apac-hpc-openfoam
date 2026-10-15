#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=26,mem=128gb
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N topo_screen_416
#PBS -j oe

# Screen topology-aware scotch multi-level domains for 416 = 4 nodes x 2 sockets x 4 NUMA x 13 cores
# (normalsr: 2 sockets, 8 NUMA domains, 13 cores each) with decomposePar -dry-run (writes nothing).
# Intended with COMPACT pinning (rank i on core i, consecutive ranks share a NUMA domain), so
# scotch's lower levels map onto sockets/NUMA domains. OpenFOAM converts domainWeights to integer
# link costs w/min(w), so powers of two keep the ratios exact.
#   cur       : domains (4 104)       weights (1 0.5)                  -- current (with scatter)
#   sock      : domains (4 2 52)      weights (1 0.5 0.25)             -- node/socket
#   numa      : domains (4 2 4 13)    weights (1 0.5 0.25 0.125)       -- node/socket/NUMA, steep
#   numaflat  : domains (4 2 4 13)    weights (1 0.5 0.5 0.5)          -- node/socket/NUMA, flat below node
# Stats show balance only (cells / neighbours / proc faces per rank), not locality.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/decomp-screen/topo416_${PBS_JOBID%%.*}
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

for v in cur sock numa numaflat; do
    cp system/include/caseDefinition.team system/include/caseDefinition
    cp system/decomposeParDict.team system/decomposeParDict
    case $v in
        cur)      d='(4 104)';    w='(1 0.5)' ;;
        sock)     d='(4 2 52)';   w='(1 0.5 0.25)' ;;
        numa)     d='(4 2 4 13)'; w='(1 0.5 0.25 0.125)' ;;
        numaflat) d='(4 2 4 13)'; w='(1 0.5 0.5 0.5)' ;;
    esac
    sed -i -e "s/^nDomains .*/nDomains $d;/" -e "s/^nDomainWeights .*/nDomainWeights $w;/" system/include/caseDefinition
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
lines = ['# topology-aware scotch domains screen, 416 ranks, occDrivAerStaticMesh 65M (decomposePar -dry-run)', '',
         '| Variant | Cells max/median | Neigh min / median / max | Faces min / median / max | Faces max/median |',
         '|---|---:|---:|---:|---:|']
for v in ('cur', 'sock', 'numa', 'numaflat'):
    s = {m[0]: tuple(map(int, m[1:])) for m in pat.findall((root / f'dryrun_{v}.log').read_text())}
    c, n, f = s['Cells'], s['Neigh'], s['Faces']
    lines.append(f'| {v} | {100*c[1]/c[2]:.2f}% | {n[0]} / {n[2]} / **{n[1]}** | {f[0]} / {f[2]} / **{f[1]}** | {100*f[1]/f[2]:.1f}% |')
lines += ['', 'cur = current 4-node setting (domains (4 104), domainWeights (1 0.5)). Lower worst-case Neigh/Faces = less MPI imbalance.',
          'Not shown: socket/NUMA locality -- the reason for the topology variants; that needs a timed run with compact pinning.']
(root / 'REPORT.md').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines))
EOF
echo "SCREEN_COMPLETE root=$ROOT"
