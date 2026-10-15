#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=40,mem=190gb
#PBS -l walltime=00:40:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N r1_decomp
#PBS -j oe

# Decomposition search round 1, step 2: real decompositions of the top-3 scored candidates
# (decomp-search/r1_416_180996035) for the paired 4-node test. Same procedure as case416_w05
# (team case, team 4-node dicts, site OpenFOAM/2512 decomposePar, collated per node, restore0Dir);
# only Basic-tier coefficients differ (caseDefinition nDomains/nDomainWeights, scotchCoeffs strategy).
#   sock211 : nDomains (4 2 52), nDomainWeights (2 1 1)
#   w5      : nDomains (4 104),  nDomainWeights (5 1)
#   sbal01  : nDomains (4 104),  nDomainWeights (1 0.5) + scotch strategy built with
#             SCOTCH_STRATBALANCE, imbalance 0.01 (tools/decompScore/bin/scotchStrat)
# The three run concurrently (~43 GB each). Each case's per-rank table is checked against the
# decompScore prediction.

set -euo pipefail
TEAM=/scratch/ya41/apac-2026/openfoam
CASE_SRC=$TEAM/run/case/occDrivAerStaticMesh
MESH=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh/constant/polyMesh
T=/scratch/ya41/kn3896/openfoam/tools/decompScore
SEARCH=/scratch/ya41/kn3896/openfoam/contest-results/decomp-search
SCORED=$SEARCH/r1_416_180996035/ranks
fh=(-fileHandler collated -ioRanks '(0 104 208 312)')
exec > >(tee "$SEARCH/r1_decompose_${PBS_JOBID%%.*}.log") 2>&1

STRAT=$("$T/bin/scotchStrat" 0x4 416 0.01)
setup() {   # $1 name, $2 nDomains, $3 nDomainWeights, $4 strategy ('' = none)
    local C=$SEARCH/case416_$1
    if [[ -e $C/system/decomposeParDict ]]; then echo "REUSING setup $C"; return 0; fi
    mkdir -p "$C/constant" "$C/logFiles"
    cp -a "$CASE_SRC/0.orig" "$CASE_SRC/system" "$C/"
    cp -a "$CASE_SRC"/constant/{transportProperties,turbulenceProperties} "$C/constant/"
    ln -s "$MESH" "$C/constant/polyMesh"
    cp "$TEAM/dicts/4node/caseDefinition"   "$C/system/include/caseDefinition"
    cp "$TEAM/dicts/4node/decomposeParDict" "$C/system/decomposeParDict"
    cp "$TEAM/dicts/4node/controlDict"      "$C/system/controlDict"
    sed -i "s/^nDomains .*/nDomains $2;/; s/^nDomainWeights .*/nDomainWeights $3;/" "$C/system/include/caseDefinition"
    grep -q "^nDomains $2;" "$C/system/include/caseDefinition"
    grep -q "^nDomainWeights $3;" "$C/system/include/caseDefinition"
    if [[ -n $4 ]]; then
        S="$4" awk '{print} /^[[:space:]]*domainWeights/ {print "    strategy \"" ENVIRON["S"] "\";"}' \
            "$C/system/decomposeParDict" > "$C/system/decomposeParDict.new"
        mv "$C/system/decomposeParDict.new" "$C/system/decomposeParDict"
        grep -qF "strategy \"$4\";" "$C/system/decomposeParDict"
    fi
}
setup sock211 '(4 2 52)' '(2 1 1)' ''
setup w5      '(4 104)'  '(5 1)'   ''
setup sbal01  '(4 104)'  '(1 0.5)' "$STRAT"

set +eu
module purge
module load OpenFOAM/2512
. "$WM_PROJECT_DIR/bin/tools/RunFunctions"   # provides restore0Dir
set -eu
type restore0Dir > /dev/null
pids=()
for n in sock211 w5 sbal01; do
    (
        cd "$SEARCH/case416_$n"
        [[ $(foamDictionary -expand system/decomposeParDict | awk '/^numberOfSubdomains/ {sub(";","",$2); print $2}') == 416 ]]
        # an earlier attempt finished sock211 and cut the others off mid-way (-force rewrites them)
        if ! grep -q '^End' logFiles/10_decomposePar.log 2>/dev/null; then
            decomposePar -constant -force "${fh[@]}" > logFiles/10_decomposePar.log 2>&1
            grep -q '^End' logFiles/10_decomposePar.log
        fi
        restore0Dir -processor > logFiles/15_restore0Dir.log 2>&1
    ) &
    pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done

# per-rank table of each real decomposition must equal the decompScore prediction
module load python3/3.12.1 > /dev/null 2>&1
declare -A SCORE=([sock211]=sock211 [w5]=w5 [sbal01]=s_bal01)
for n in sock211 w5 sbal01; do
    python3 -I - "$SEARCH/case416_$n/logFiles/10_decomposePar.log" "$SCORED/${SCORE[$n]}.ranks" "$n" <<'PY'
import re, sys
log = open(sys.argv[1]).read()
dp = [(int(p), int(c), int(f), int(b), int(n)) for p, c, n, f, b in re.findall(
    r'^Processor (\d+)\n\s+Number of cells = (\d+)\n.*?Number of processor patches = (\d+)\n\s+'
    r'Number of processor faces = (\d+)\n\s+Number of boundary faces = (\d+)', log, re.S | re.M)]
pred = [tuple(int(x) for x in l.split()[:5]) for l in open(sys.argv[2]) if not l.startswith('#')]
print(f'CHECK {sys.argv[3]}: ranks {len(dp)} matches_decompScore={"yes" if dp == pred else "NO"}')
PY
done
echo "R1_DECOMPOSE_SUCCESS"
