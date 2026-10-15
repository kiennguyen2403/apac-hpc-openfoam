#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=40,mem=190gb,jobfs=2gb
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N decomp_search
#PBS -j oe

# Decomposition search, round 1 (4 nodes = 416 ranks, occDrivAerStaticMesh 65M, Basic-tier
# method coefficients only). decompScore (tools/decompScore, one forked child per candidate, site OpenFOAM/2512 + scotch 6.1.0,
# the build that decomposed case416_w05) reads the mesh once per instance and writes per-rank
# cells / processor faces / boundary faces / neighbours for each candidate; nothing is
# decomposed to disk. score_decomp.py ranks candidates by the predicted slowest rank
# (per-rank model fitted on IPM job 180805820).
# Candidates:
#   wK        scotch (4 104), link costs (K 1)  [OpenFOAM truncates domainWeights to integer
#             ratios: (1 0.5) == (2 1); w2 = current, w100 = old (1 0.01)]
#   s_*       scotch (4 104)/(2 1) with a scotch-built strategy string (flags, imbalance)
#   plain*    scotch without domains
#   sock*     scotch (4 2 52) with link costs
#   h_AxBxC   hierarchical n (A B C), order xyz, delta 0.001
# Self-checks: w2 must reproduce case416_w05 rank by rank; w100 vs w2 is the model check
# (measured: w2 is 1.15% faster, job 180782454).

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
NINST=6          # only used to order the candidate list
NPAR=6           # candidates decomposed concurrently (forked children share the mesh)
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh
T=/scratch/ya41/kn3896/openfoam/tools/decompScore
IPM=/scratch/ya41/kn3896/openfoam/contest-results/ipm-site/ipm_4n_180805820/kn3896.1791451559.186205.ipm.xml
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/decomp-search/r1_416_${PBS_JOBID%%.*}
W=$PBS_JOBFS/work

mkdir -p "$ROOT/ranks" "$ROOT/logs" "$W/case/constant" "$W/case/system" "$W/dicts"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "SEARCH_FAILED line=$LINENO status=$?"' ERR
[[ $(grep -a -m1 -o 'nCells:[0-9]*' "$MESH/constant/polyMesh/owner") == nCells:65334765 ]]
ln -s "$MESH/constant/polyMesh" "$W/case/constant/polyMesh"
cat > "$W/case/system/controlDict" <<'DICT'
FoamFile { version 2.0; format ascii; class dictionary; object controlDict; }
application none; startFrom startTime; startTime 0; stopAt endTime; endTime 1; deltaT 1;
writeControl timeStep; writeInterval 1;
DICT

# --- candidates ------------------------------------------------------------------------
strat() { "$T/bin/scotchStrat" "$1" 416 "$2"; }
python3 -I - "$W/dicts" "$NINST" \
    "$(strat 0x0 0.01)" "$(strat 0x0 0.05)" "$(strat 0x1 0.01)" "$(strat 0x1 0.05)" \
    "$(strat 0x4 0.01)" "$(strat 0x5 0.01)" "$(strat 0x100 0.01)" <<'PY'
import sys
from pathlib import Path
d, ninst = Path(sys.argv[1]), int(sys.argv[2])
s_def01, s_def05, s_qual01, s_qual05, s_bal01, s_qualbal01, s_rec01 = sys.argv[3:10]
head = 'FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }\nnumberOfSubdomains 416;\n'
cands = []   # (name, body, cost_s, priority)  priority: 0 first, 2 last (strategy strings)
for k in (2, 100, 1, 3, 4, 5, 7, 10):
    cands.append((f'w{k}', f'method scotch;\ncoeffs {{ domains (4 104); domainWeights ({k} 1); }}\n', 270, 0 if k in (2, 100) else 1))
for name, s in [('s_def01', s_def01), ('s_def05', s_def05), ('s_qual01', s_qual01), ('s_qual05', s_qual05),
                ('s_bal01', s_bal01), ('s_qualbal01', s_qualbal01), ('s_rec01', s_rec01)]:
    cands.append((name, f'method scotch;\ncoeffs {{ domains (4 104); domainWeights (2 1); strategy "{s}"; }}\n', 270, 2))
cands.append(('plain', 'method scotch;\ncoeffs { }\n', 270, 1))
cands.append(('plain_qual01', f'method scotch;\ncoeffs {{ strategy "{s_qual01}"; }}\n', 270, 2))
cands.append(('sock421', 'method scotch;\ncoeffs { domains (4 2 52); domainWeights (4 2 1); }\n', 270, 1))
cands.append(('sock211', 'method scotch;\ncoeffs { domains (4 2 52); domainWeights (2 1 1); }\n', 270, 1))
divs = [x for x in range(1, 417) if 416 % x == 0]
for a in divs:
    for b in divs:
        if (416 // a) % b == 0:
            c = 416 // (a * b)
            cands.append((f'h_{a}x{b}x{c}', f'method hierarchical;\ncoeffs {{ n ({a} {b} {c}); order xyz; delta 0.001; }}\n', 30, 1))
# longest-processing-time assignment; within an instance: priority, then cost (long first)
load, lists = [0] * ninst, [[] for _ in range(ninst)]
for c in sorted(cands, key=lambda c: (c[3] == 0, c[2]), reverse=True):
    i = load.index(min(load)); load[i] += c[2]; lists[i].append(c)
for i, lst in enumerate(lists):
    lst.sort(key=lambda c: (c[3], -c[2]))
    lines = []
    for name, body, _, _ in lst:
        (d / f'{name}').write_text(head + body)
        lines.append(f'{name} {d / name}')
    (d / f'inst{i}.txt').write_text('\n'.join(lines) + '\n')
print(f'CANDIDATES {len(cands)} instances {ninst} est_load_s {load}')
PY
for f in "$W"/dicts/inst*.txt; do while read -r n p; do printf '### %s\n' "$n"; cat "$p"; done < "$f"; done > "$ROOT/candidates.txt"

# --- run decompScore instances (site OpenFOAM/2512) ------------------------------------
set +eu
module purge
module load openmpi/4.1.7
source /apps/OpenFOAM/2512/etc/bashrc > /dev/null 2>&1
set -eu
[[ $WM_OPTIONS == linux64GccDPInt32Opt ]]
export FOAM_SIGFPE=false OMP_NUM_THREADS=1
echo "SEARCH_START utc=$(date -u +%FT%TZ)"
cat "$W"/dicts/inst*.txt > "$W/dicts/all.txt"
fails=0
"$T/bin/decompScore" -case "$W/case" -candidates "$W/dicts/all.txt" -outDir "$ROOT/ranks" -nParallel "$NPAR" \
    > "$ROOT/logs/decompScore.log" 2>&1 || fails=1
grep -h '^CANDIDATE_FAILED' "$ROOT/logs/decompScore.log" || true
echo "SEARCH_END utc=$(date -u +%FT%TZ) decompScore_failed=$fails"
grep -h '^CANDIDATE ' "$ROOT"/logs/decompScore.log | sort -k2,2 > "$ROOT/candidates_summary.txt" || true
missing=$(cat "$W"/dicts/inst*.txt | while read -r n p; do [[ -f $ROOT/ranks/$n.ranks ]] || echo "$n"; done | tr '\n' ' ')
echo "MISSING_CANDIDATES: ${missing:-none}"

# --- self-checks and scoring -----------------------------------------------------------
if cmp -s <(grep -v '^#' "$ROOT/ranks/w2.ranks") <(grep -v '^#' "$T/train_case416_w05.ranks"); then
    echo "CHECK_W2_REPRODUCES_CASE416_W05: yes"
else
    echo "CHECK_W2_REPRODUCES_CASE416_W05: NO"
fi
if cmp -s <(grep -v '^#' "$ROOT/ranks/w2.ranks") <(grep -v '^#' "$ROOT/ranks/s_def01.ranks" 2>/dev/null); then
    echo "CHECK_SCOTCH_DEFAULT_EQUALS_DEFAULT_0.01: yes"
else
    echo "CHECK_SCOTCH_DEFAULT_EQUALS_DEFAULT_0.01: no"
fi
module load python3/3.12.1 > /dev/null 2>&1
python3 -I "$T/score_decomp.py" "$IPM" "$T/train_case416_w05.ranks" "$ROOT/ranks" w2 "$ROOT/REPORT.md"
echo "SEARCH_COMPLETE root=$ROOT"
