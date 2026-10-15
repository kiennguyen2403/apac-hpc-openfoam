#!/bin/bash -l

set -eEuo pipefail
trap 'rc=$?; echo "ERROR line $LINENO: $BASH_COMMAND (status $rc)" >&2' ERR

: "${PBS_JOBID:?Run this script inside a PBS job}"
: "${PBS_NODEFILE:?PBS_NODEFILE is unavailable}"
: "${PROFILE_MODE:=contest}"
: "${PROFILE_VARIANTS:=intelmpi openmpi mpich}"
: "${RESULT_ROOT:=/scratch/ya41/kn3896/openfoam/contest-results/communication}"
: "${CONTEST_CASE:=/scratch/ya41/kn3896/cases/openfoam-hpc/incompressible/simpleFoam/occDrivAerStaticMesh}"
export CONTEST_CASE

RUNNER=/home/565/kn3896/drivaer-contest-profile/run-contest-case.sh
LOADER=/home/565/kn3896/drivaer-contest-profile/load-communication-variant.sh
CUSTOM_FOAM=/scratch/ya41/kn3896/openfoam/OpenFOAM-v2512

run_smoke()
(
    local variant=$1
    local ranks=${SMOKE_RANKS:-2}
    local smoke_nodes=${SMOKE_NODES:-1}
    local case_dir=${SMOKE_CASE_ROOT:-$PBS_JOBFS}/communication-smoke-$variant
    [[ $(sort -u "$PBS_NODEFILE" | wc -l) -eq $smoke_nodes ]]
    if ((smoke_nodes > 1)); then
        [[ -n ${SMOKE_CASE_ROOT:-} && $ranks -eq $smoke_nodes ]] || {
            echo "ERROR: multi-node smoke requires shared SMOKE_CASE_ROOT and one rank per node" >&2
            exit 2
        }
    fi
    local result_dir=$RESULT_ROOT/smoke_${PBS_JOBID}/$variant

    export FOAM_VARIANT=$variant
    set +eu
    saved_err_trap=$(trap -p ERR || true)
    trap - ERR
    source "$LOADER"
    loader_status=$?
    eval "$saved_err_trap"
    set -eu
    [[ $loader_status -eq 0 ]]

    for command_name in blockMesh decomposePar simpleFoam; do
        command -v "$command_name" >/dev/null
    done

    mkdir -p "$result_dir"
    sort -u "$PBS_NODEFILE" > "$result_dir/hosts"
    if ((smoke_nodes > 1)); then
        case "$MPI_FLAVOUR" in
            openmpi) MPI_PLACEMENT=(--map-by ppr:1:node --bind-to core) ;;
            intelmpi|mpich) MPI_PLACEMENT+=(-ppn 1 -f "$result_dir/hosts") ;;
        esac
    fi
    {
        echo "variant=$variant"
        echo "ranks=$ranks"
        echo "wm_project_version=$WM_PROJECT_VERSION"
        echo "wm_options=$WM_OPTIONS"
        echo "foam_mpi=$FOAM_MPI"
        echo "mpi_arch_path=$MPI_ARCH_PATH"
        echo "simpleFoam=$(command -v simpleFoam)"
        printf 'launcher='
        printf '%q ' "${MPI_LAUNCHER[@]}" "$ranks" "${MPI_PLACEMENT[@]}"
        printf '\n'
    } > "$result_dir/metadata.txt"
    mpi=("${MPI_LAUNCHER[@]}" "$ranks" "${MPI_PLACEMENT[@]}")
    "${mpi[@]}" /bin/bash /home/565/kn3896/drivaer-contest-profile/check-communication-runtime.sh \
        simpleFoam decomposePar renumberMesh potentialFoam applyBoundaryLayer \
        > "$result_dir/runtime-check.log" 2>&1
    [[ $(awk '/^RUNTIME_SUCCESS/ {for(i=1;i<=NF;i++) if($i ~ /^host=/) print $i}' \
        "$result_dir/runtime-check.log" | sort -u | wc -l) -eq $smoke_nodes ]]
    cp -a "$CUSTOM_FOAM/tutorials/incompressible/simpleFoam/rotatingCylinders" "$case_dir"
    cp "$WM_PROJECT_DIR/etc/caseDicts/annotated/decomposeParDict" "$case_dir/system/decomposeParDict"
    foamDictionary -entry numberOfSubdomains -set "$ranks" \
        "$case_dir/system/decomposeParDict" >/dev/null
    foamDictionary -entry method -set scotch \
        "$case_dir/system/decomposeParDict" >/dev/null
    foamDictionary -entry endTime -set 2 "$case_dir/system/controlDict" >/dev/null

    blockMesh -case "$case_dir" > "$result_dir/blockMesh.log" 2>&1
    decomposePar -case "$case_dir" > "$result_dir/decomposePar.log" 2>&1
    mpi=("${MPI_LAUNCHER[@]}" "$ranks" "${MPI_PLACEMENT[@]}")
    "${mpi[@]}" simpleFoam -parallel -case "$case_dir" \
        > "$result_dir/simpleFoam.log" 2>&1

    grep -q '^End$' "$result_dir/simpleFoam.log"
    ! grep -Eq 'FOAM FATAL|Segmentation fault|MPI_ABORT' "$result_dir/simpleFoam.log"
    echo "smoke_status=passed" >> "$result_dir/metadata.txt"
    rm -rf "$case_dir"
    echo "SMOKE_SUCCESS variant=$variant"
)

case "$PROFILE_MODE" in
    smoke)
        : "${PBS_JOBFS:?Smoke test requires PBS jobfs}"
        for variant in $PROFILE_VARIANTS; do
            run_smoke "$variant"
        done
        ;;
    ipm-smoke)
        source /home/565/kn3896/drivaer-contest-profile/ipm-functions.sh
        run_ipm_smoke
        ;;
    pilot)
        source /home/565/kn3896/drivaer-contest-profile/pilot-functions.sh
        run_pilot
        ;;
    contest)
        : "${NODES:?Set NODES=1 or NODES=2}"
        case "$NODES" in
            1) expected_end=40 ;;
            2) expected_end=80 ;;
            *) echo "ERROR: communication profiling supports NODES=1 or 2" >&2; exit 2 ;;
        esac
        : "${END_TIME:=$expected_end}"
        [[ $END_TIME -eq $expected_end ]]
        for variant in $PROFILE_VARIANTS; do
            export FOAM_VARIANT=$variant
            export RESULT_LABEL=${NODES}node_${END_TIME}iter_${variant}_${PBS_JOBID}
            export KEEP_CASE=0 PROFILE_IPM=0
            "$RUNNER"
        done
        ;;
    *)
        echo "ERROR: PROFILE_MODE must be smoke, ipm-smoke, pilot, or contest" >&2
        exit 2
        ;;
esac
