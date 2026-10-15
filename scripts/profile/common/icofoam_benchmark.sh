#!/usr/bin/env bash

ROOT=/home/565/kn3896/scratch/openfoam
FOAM_DIR="$ROOT/OpenFOAM-dev"
BASE_CASE="$FOAM_DIR/tutorials/incompressible/icoFoam/cavity/cavity"
NPROCS="${PBS_NCPUS:-8}"
MESH_N="${MESH_N:-256}"
REPEATS="${REPEATS:-3}"
END_TIME="${END_TIME:-0.005}"

bench_init()
{
    RUN_DIR="$ROOT/profiles/$1_${PBS_JOBID:-manual}"
    LOG_DIR="$ROOT/logs/profiles/$1"
    mkdir -p "$RUN_DIR" "$LOG_DIR"
    exec > >(tee "$LOG_DIR/${PBS_JOBID:-manual}.log") 2>&1
    trap 'rc=$?; echo "ERROR: line $LINENO: $BASH_COMMAND (exit $rc)" >&2' ERR
    [[ "$NPROCS" == 8 && -d "$BASE_CASE" ]]
}

bench_run()
{
    local label=$1 mpi_exec=$2 case_dir="$RUN_DIR/$1"
    for cmd in blockMesh decomposePar foamDictionary icoFoam; do
        command -v "$cmd" >/dev/null || { echo "ERROR: $label missing $cmd"; return 1; }
    done
    cp -a "$BASE_CASE" "$case_dir"
    sed -i "s/(20 20 1)/($MESH_N $MESH_N 1)/" "$case_dir/system/blockMeshDict"
    foamDictionary -entry numberOfSubdomains -set "$NPROCS" "$case_dir/system/decomposeParDict" >/dev/null
    foamDictionary -entry method -set hierarchical "$case_dir/system/decomposeParDict" >/dev/null
    foamDictionary -entry coeffs -remove "$case_dir/system/decomposeParDict" >/dev/null
    foamDictionary -entry hierarchicalCoeffs -set '{ n (4 2 1); order xyz; }' "$case_dir/system/decomposeParDict" >/dev/null
    foamDictionary -entry deltaT -set 0.00001 "$case_dir/system/controlDict" >/dev/null
    foamDictionary -entry endTime -set "$END_TIME" "$case_dir/system/controlDict" >/dev/null
    foamDictionary -entry writeInterval -set 500 "$case_dir/system/controlDict" >/dev/null
    blockMesh -case "$case_dir" > "$case_dir/log.blockMesh" 2>&1
    decomposePar -case "$case_dir" > "$case_dir/log.decomposePar" 2>&1
    {
        echo "label=$label"
        echo "version=${WM_PROJECT_VERSION:-unknown}"
        echo "platform=${WM_OPTIONS:-unknown}"
        echo "mpi=${FOAM_MPI:-${WM_MPLIB:-unknown}}"
        "${WM_CXX:-c++}" --version 2>&1 | head -1 || true
        "$mpi_exec" --version 2>&1 | head -2 || true
    } > "$case_dir/build-info.txt"
    for ((run=1; run<=REPEATS; run++)); do
        echo "$label: run $run/$REPEATS"
        /usr/bin/time -f 'elapsed_seconds=%e\nuser_seconds=%U\nsystem_seconds=%S\nmax_rss_kb=%M' \
            -o "$case_dir/time.$run.txt" \
            "$mpi_exec" -n "$NPROCS" icoFoam -parallel -case "$case_dir" \
            > "$case_dir/log.icoFoam.$run" 2>&1
        grep -q '^End$' "$case_dir/log.icoFoam.$run"
        ! grep -Eq 'FOAM FATAL|Segmentation fault' "$case_dir/log.icoFoam.$run"
    done
}

bench_summary()
{
    local label value fastest row
    local -a rows=()
    for label; do
        value=$(awk -F= '/^elapsed_seconds=/{s+=$2;n++} END{if(n)printf "%.6f",s/n}' "$RUN_DIR/$label"/time.*.txt)
        rows+=("$label:$value")
    done
    fastest=$(printf '%s\n' "${rows[@]#*:}" | sort -n | head -1)
    printf 'variant\taverage_seconds\trelative_to_fastest\n' > "$RUN_DIR/summary.tsv"
    for row in "${rows[@]}"; do
        label=${row%%:*}; value=${row#*:}
        printf '%s\t%s\t%s\n' "$label" "$value" \
            "$(awk -v v="$value" -v f="$fastest" 'BEGIN{printf "%.4f",v/f}')" \
            >> "$RUN_DIR/summary.tsv"
    done
    column -t -s $'\t' "$RUN_DIR/summary.tsv" || cat "$RUN_DIR/summary.tsv"
    echo "Results: $RUN_DIR"
}
