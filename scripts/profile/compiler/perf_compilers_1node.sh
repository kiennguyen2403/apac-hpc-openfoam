#!/bin/bash
#PBS -P ya41
#PBS -q normalsr
#PBS -l ncpus=104,mem=400gb
#PBS -l walltime=01:00:00
#PBS -l storage=scratch/ya41
#PBS -l wd
#PBS -N perf_cc_1n
#PBS -j oe

# perf profile of vanilla / GCC / Clang / ICX on occDrivAerStaticMesh 65M, 1 node, 104 ranks, 40 steps,
# all with the SAME MPI library and runtime settings: Open MPI 4.1.7 (the version module
# OpenFOAM/2512 is built with), team scatter CPU order pinned via --cpu-list, default collectives.
#   vanilla      : module OpenFOAM/2512 (GCC 8.5, site build)
#   gcc/clang/icx: custom v2512 OptSPROmpi builds (-O3 + Sapphire Rapids target, -fp-model precise
#                  for clang/icx), built against openmpi/4.1.7
# Reuses the prepared 1-node scotch case (no decomposition).
# Per compiler: 3-step warmup (checks env, builds the coded timer for that platform), then one
# 40-step perf run: perf stat on all ranks except 0 26 52 78, perf record (cycles:u, 499 Hz) on those.

set -euo pipefail
: "${PBS_JOBID:?submit through PBS}"
NP=104
FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512
TEAM=/scratch/ya41/apac-2026/openfoam
CASE=/scratch/ya41/kn3896/openfoam/contest-results/vtune/drivaer_1node_icx_180728663/case
ROOT=/scratch/ya41/kn3896/openfoam/contest-results/perf-compilers/drivaer_1node_${PBS_JOBID%%.*}
fh="-fileHandler hostCollated"
RECORD_RANKS="0 26 52 78"
STEPS=40
COMPILERS="${PERF_COMPILERS:-vanilla gcc clang icx}"

mkdir -p "$ROOT"
exec > >(tee "$ROOT/job.log") 2>&1
trap 'echo "PERFCC_FAILED line=$LINENO status=$?"' ERR
cd "$ROOT"
cp "$TEAM/run/1node/intel-scatter-cpus.awk" .
for l in 10_decomposePar 20_renumberMesh 30_potentialFoam 40_applyBoundaryLayer; do grep -q '^End' "$CASE/../logFiles/$l.log"; done
cmp "$CASE/system/fvSolution" "$CASE/system/fvSolution.fixedIter"
if [[ -d $CASE/postProcessing ]]; then mv "$CASE/postProcessing" "$ROOT/postProcessing.before"; fi
cat /proc/sys/kernel/perf_event_paranoid > perf_event_paranoid.txt
lscpu > lscpu.txt

# Per-rank pinning via taskset: Gadi's Open MPI config forces a mapping policy that conflicts
# with rankfiles, and --cpu-list ... :ordered does not keep the list order.
cat > pin.sh <<'EOF'
#!/bin/bash
# Pin this rank to the LOCAL_RANK-th CPU of PIN_LIST (team scatter order), then exec the command.
IFS=, read -ra cpus <<< "${PIN_LIST:?}"
exec taskset -c "${cpus[${OMPI_COMM_WORLD_LOCAL_RANK:?}]}" "$@"
EOF
chmod +x pin.sh

cat > perf-wrap.sh <<'EOF'
#!/bin/bash
r=${OMPI_COMM_WORLD_RANK:?no MPI rank variable}
case " $RECORD_RANKS " in
    *" $r "*) exec perf record -q -F 499 -e cycles:u -o "$PERF_DIR/record.$r.data" -- "$@" ;;
    *) exec perf stat -x, -o "$PERF_DIR/stat.$r.csv" \
           -e cycles:u,instructions:u,branch-misses:u,cache-misses:u \
           -e cpu/event=0xc7,umask=0x01,name=fp_scalar_double/u \
           -e cpu/event=0xc7,umask=0x04,name=fp_128b_packed_double/u \
           -e cpu/event=0xc7,umask=0x10,name=fp_256b_packed_double/u \
           -e cpu/event=0xc7,umask=0x40,name=fp_512b_packed_double/u \
           -- "$@" ;;
esac
EOF
chmod +x perf-wrap.sh

load_variant() {   # $1 = vanilla|gcc|clang|icx ; runs inside a subshell
    local comp cc cxx mods flags
    set +eu
    module purge
    case $1 in
        vanilla)
            module load OpenFOAM/2512
            export WM_COMPILER_TYPE=system WM_CC=gcc WM_CXX=g++ CC=gcc CXX=g++
            ;;
        gcc|clang|icx|icxfma)
            local opt=OptSPROmpi
            case $1 in
                gcc)   mods="gcc/14.2.0 binutils/2.44"; comp=Gcc;   cc=gcc;   cxx=g++;     flags="-march=sapphirerapids" ;;
                clang) mods="clang/19.1.6";             comp=Clang; cc=clang; cxx=clang++; flags="-march=sapphirerapids -ffp-model=precise" ;;
                icx)   mods="intel-compiler-llvm/2025.3.2"; comp=Icx; cc=icx;  cxx=icpx;    flags="-xSAPPHIRERAPIDS -fp-model=precise" ;;
                icxfma) mods="intel-compiler-llvm/2025.3.2"; comp=Icx; cc=icx; cxx=icpx;    flags="-xSAPPHIRERAPIDS -fp-model=precise -ffp-contract=fast"; opt=OptSPRfmaOmpi ;;
            esac
            module load $mods openmpi/4.1.7 metis/5.1.0 fftw3/3.3.10 boost/1.84.0
            export CC=$cc CXX=$cxx OMPI_CC=$cc OMPI_CXX=$cxx FOAM_EXTRA_CFLAGS="$flags" FOAM_EXTRA_CXXFLAGS="$flags"
            source "$FOAM/etc/bashrc" WM_COMPILER_TYPE=system WM_COMPILER=$comp WM_CC=$cc WM_CXX=$cxx \
                WM_MPLIB=SYSTEMOPENMPI WM_COMPILE_OPTION=$opt
            export WM_THIRD_PARTY_DIR=/apps/OpenFOAM/2512/ThirdParty-v2512
            local pfx; pfx=$(ompi_info --path prefix --parsable); export MPI_ARCH_PATH=${pfx#path:prefix:} FOAM_MPI=sys-openmpi
            ;;
    esac
    set -eu
    case $1 in
        vanilla) [[ $WM_OPTIONS == linux64GccDPInt32Opt ]] || { echo "ERROR: WM_OPTIONS=$WM_OPTIONS"; return 1; } ;;
        *)       [[ $WM_OPTIONS == linux64${comp}DPInt32$opt ]] || { echo "ERROR: WM_OPTIONS=$WM_OPTIONS"; return 1; } ;;
    esac
    [[ $(mpirun --version 2>&1 | head -1) == *"Open MPI) 4.1.7"* ]] || { echo "ERROR: not Open MPI 4.1.7"; return 1; }
    export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 FOAM_SIGFPE=true
    # Team scatter CPU order, rank i pinned to the i-th CPU of the list (same placement as the Intel MPI runs)
    SCATTER=$(lscpu -p=CPU,CORE,SOCKET,NODE,ONLINE | sort -t, -k3,3n -k4,4n -k1,1n | gawk -f "$ROOT/intel-scatter-cpus.awk")
    export PIN_LIST=$SCATTER PERF_DIR=${PERF_DIR:-none} RECORD_RANKS   # -x needs them set
    MPI_RUN=(mpirun -np $NP --bind-to none -x PIN_LIST -x PERF_DIR -x RECORD_RANKS "$ROOT/pin.sh")
}

run() {   # $1 compiler, $2 warmup|perf ; failures are recorded and the job continues
    local c=$1 label=$2 dir=$ROOT/run_${1}_$2
    mkdir "$dir"
    (
        load_variant "$c" || exit 2
        local steps=$STEPS wrap=() status=0
        if [[ $label == warmup ]]; then
            steps=3
            {
                echo "compiler=$c wm_options=$WM_OPTIONS"
                echo "simpleFoam=$(command -v simpleFoam)"
                echo "compiler_version=$($CXX --version | head -1)"
                echo "mpi=$(mpirun --version 2>&1 | head -1)"
                echo "flags=$(wmake -show-cxxflags)"
            } > "$dir/build-metadata.txt"
            ldd "$(command -v simpleFoam)" > "$dir/simpleFoam-ldd.txt" || exit 3
            if grep -q 'not found' "$dir/simpleFoam-ldd.txt"; then exit 3; fi
            # /apps/openmpi/4.1.7 is a link to /apps/openmpi-mofed5.8-pbs2021.1/4.1.7 (same install)
            grep -qE 'openmpi(-mofed5[.]8-pbs2021[.]1)?/4[.]1[.]7/lib/libmpi' "$dir/simpleFoam-ldd.txt" || { echo "ERROR: not linked to Open MPI 4.1.7"; exit 3; }
            # placement check: every rank must sit on its scatter-list CPU
            "${MPI_RUN[@]}" bash -c 'echo "$OMPI_COMM_WORLD_RANK $(taskset -cp $$ | awk "{print \$NF}")"' | sort -n > "$dir/placement.txt"
            [[ $(awk '{print $2}' "$dir/placement.txt" | paste -sd,) == "$PIN_LIST" ]] || { echo "ERROR: placement differs from scatter list"; exit 3; }
            # coded timer: build for this platform if missing (OpenFOAM only rebuilds on source change)
            local td=$CASE/dynamicCode/wallClockTimeStatistics
            if ! ls "$td"/../platforms/$WM_OPTIONS/lib/libwallClockTimeStatistics_*.so > /dev/null 2>&1; then
                (cd "$td" && wmake libso) > "$dir/timer-build.log" 2>&1 || exit 3
            fi
        else
            export PERF_DIR=$dir RECORD_RANKS
            wrap=("$ROOT/perf-wrap.sh")
        fi
        cd "$CASE"
        cp "$TEAM/dicts/1node/controlDict" system/controlDict
        sed -i "s/^endTime .*/endTime           $steps;/" system/controlDict
        [[ $(sed -n 's/^endTime *\([0-9]*\);.*/\1/p' system/controlDict) == $steps ]] || exit 4
        echo "RUN_START compiler=$c run=$label utc=$(date -u +%FT%TZ)"
        timeout --kill-after=30s 1500s "${MPI_RUN[@]}" "${wrap[@]}" simpleFoam -parallel $fh > "$dir/simpleFoam.log" 2>&1 || status=$?
        cp "$TEAM/dicts/1node/controlDict" system/controlDict
        if [[ -d postProcessing ]]; then mv postProcessing "$dir/"; fi
        [[ $status == 0 ]] || { echo "RUN_FAILED compiler=$c run=$label status=$status"; exit 5; }
        grep -qx End "$dir/simpleFoam.log" || exit 6
        grep -q 'Average wall-clock time per time step' "$dir/simpleFoam.log" || { echo "RUN_FAILED compiler=$c run=$label: no timer"; exit 7; }
        if [[ $label == perf ]]; then
            for r in $RECORD_RANKS; do
                perf report -i "$dir/record.$r.data" --stdio --sort dso,sym --percent-limit 0.3 > "$dir/report.$r.txt" 2>/dev/null || true
            done
        fi
        local cd=$(awk '!/^#/ {last=$2} END {print last}' "$dir"/postProcessing/forceCoeffsAll/*/coefficient*.dat)
        echo "RUN_SUCCESS compiler=$c run=$label Cd_step$steps=$cd $(grep 'Average wall-clock time per time step' "$dir/simpleFoam.log" | tail -1)"
    )
}

OK=()
for c in $COMPILERS; do
    case $c in
        vanilla) bin=/apps/OpenFOAM/2512/platforms/linux64GccDPInt32Opt/bin/simpleFoam ;;
        gcc)     bin=$FOAM/platforms/linux64GccDPInt32OptSPROmpi/bin/simpleFoam ;;
        clang)   bin=$FOAM/platforms/linux64ClangDPInt32OptSPROmpi/bin/simpleFoam ;;
        icx)     bin=$FOAM/platforms/linux64IcxDPInt32OptSPROmpi/bin/simpleFoam ;;
        icxfma)  bin=$FOAM/platforms/linux64IcxDPInt32OptSPRfmaOmpi/bin/simpleFoam ;;
    esac
    if [[ ! -x $bin ]]; then echo "SKIP compiler=$c (build missing: $bin)"; continue; fi
    if run "$c" warmup; then OK+=("$c"); else echo "COMPILER_DROPPED compiler=$c (warmup failed)"; fi
done
for c in "${OK[@]}"; do run "$c" perf || echo "RUN_FAILED_RECORDED compiler=$c"; done
cp "$TEAM/dicts/1node/controlDict" "$CASE/system/controlDict"

python3 - "$ROOT" "$RECORD_RANKS" "${OK[@]}" <<'EOF'
import csv, re, sys
from pathlib import Path
root, rec, comps = Path(sys.argv[1]), sys.argv[2].split(), sys.argv[3:]
pat = re.compile(r'Average wall-clock time per time step\s*=\s*([0-9.eE+-]+)')
W = {'fp_scalar_double': 1, 'fp_128b_packed_double': 2, 'fp_256b_packed_double': 4, 'fp_512b_packed_double': 8}
def stat(d):
    tot, n = {}, 0
    for p in sorted(d.glob('stat.*.csv')):
        n += 1
        for row in csv.reader(p.open()):
            if len(row) < 3 or not row[0].strip() or row[0].startswith('#'): continue
            try: v = float(row[0])
            except ValueError: continue
            k = row[2].split('/')[0].replace(':u', '')
            tot[k] = tot.get(k, 0) + v
    return n, tot
def hot(f, k=12):
    out = []
    for line in (f.read_text().splitlines() if f.exists() else []):
        m = re.match(r'\s+([0-9.]+)%\s+(\S+)\s+\[\.\]\s+(.*)$', line)
        if m: out.append((float(m[1]), m[2], m[3].strip()))
    return out[:k]
L = ['# perf: vanilla vs GCC vs Clang vs ICX, same Open MPI 4.1.7 + settings, occDrivAerStaticMesh 65M, 1 node x 104, 40 steps', '',
     'vanilla = module OpenFOAM/2512 (GCC 8.5); others = custom v2512 OptSPROmpi builds; all Open MPI 4.1.7, team scatter CPU order. Times are from the',
     'perf run itself (perf stat on 100 ranks + perf record on 4 ranks; overhead ~0.5%, equal for all builds).', '',
     '| Compiler | s/step (perf run) | Cd step 40 | Instructions (1e12) | IPC | FP elements (1e12) | 512b | 256b | 128b | scalar | cache-misses (1e9) |',
     '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for c in comps:
    d = root / f'run_{c}_perf'
    log = d / 'simpleFoam.log'
    t = pat.findall(log.read_text())[-1] if log.exists() and pat.findall(log.read_text()) else 'n/a'
    rows = [l.split() for p in (d / 'postProcessing/forceCoeffsAll').glob('*/coefficient*.dat') for l in p.read_text().splitlines() if l.strip() and not l.startswith('#')]
    cd = rows[-1][1] if rows else 'n/a'
    n, s = stat(d)
    if 'instructions' in s and 'cycles' in s:
        e = {k: s.get(k, 0) * w for k, w in W.items()}; T = sum(e.values()) or float('nan')
        L.append(f"| {c} | {t} | {float(cd):.5f} | {s['instructions']/1e12:.3f} | {s['instructions']/s['cycles']:.2f} | {T/1e12:.3f} | "
                 f"{100*e['fp_512b_packed_double']/T:.1f}% | {100*e['fp_256b_packed_double']/T:.1f}% | {100*e['fp_128b_packed_double']/T:.1f}% | "
                 f"{100*e['fp_scalar_double']/T:.1f}% | {s.get('cache-misses',0)/1e9:.1f} |")
    else:
        L.append(f'| {c} | {t} | {cd} | perf stat missing | | | | | | | |')
for c in comps:
    L += ['', f'## Hotspots {c} (perf record cycles:u, rank {rec[0]})', '', '| % cycles | Library | Symbol |', '|---:|---|---|']
    L += [f'| {p:.1f} | {lib} | `{sym[:90]}` |' for p, lib, sym in hot(root / f'run_{c}_perf' / f'report.{rec[0]}.txt')] or ['| | n/a | |']
(root / 'REPORT.md').write_text('\n'.join(L) + '\n')
print('\n'.join(L[:12]))
EOF
echo "PERFCC_COMPLETE root=$ROOT"
