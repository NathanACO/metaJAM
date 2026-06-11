#!/bin/bash
# Usage: ./kraken2_parallele_memory_mapping.sh <list_fq> <outdir> <k2_db> <label> [old_label] [threads] [parallel_jobs]

# --- Argument parsing ---
list_fq="${1:?Error: list_fq not provided}"
outdir="${2:?Error: outdir not provided}"
k2_db="${3:?Error: k2_db not provided}"
label="${4:?Error: label not provided}"
old_label="${5:-}"          # optional, default empty
threads="${6:-1}"           # optional, default 1
parallel_jobs="${7:-4}"     # optional, default 4

export outdir k2_db label old_label threads

# --- Validate inputs ---
if [[ ! -r "$list_fq" ]]; then
    echo "Error: $list_fq not found or not readable."
    exit 1
fi
if [[ ! -d "$k2_db" ]]; then
    echo "Error: database path $k2_db not found."
    exit 1
fi

# --- Function ---
run_kraken2() {
    local reads="$1"
    local data_from
    data_from=$(dirname "$reads")
    local name
    name=$(basename "$reads" \
    | sed 's/\.fastq\.gz$//' \
    | sed "s/_${old_label}$//" \
    | sed 's/_unclas$//' \
    | sed 's/\.merged$//')

    local tmp_db="/tmp/$(basename "$k2_db")"

    kraken2 \
        --db "$tmp_db" \
        --report "$outdir/${name}_${label}_report.txt" \
        --report-minimizer-data \
        --gzip-compressed \
        --threads "$threads" \
        --output "$outdir/${name}_${label}_output.txt" \
        --classified-out "$outdir/${name}_${label}_clas.fastq" \
        --unclassified-out "$outdir/${name}_${label}_unclas.fastq" \
        --memory-mapping \
        "$reads"

    pigz "$outdir/${name}_${label}_clas.fastq" && \
    pigz "$outdir/${name}_${label}_unclas.fastq" && \
    pigz "$outdir/${name}_${label}_output.txt" && \
    pigz "$outdir/${name}_${label}_report.txt"
}
export -f run_kraken2

# --- Copy DB to /tmp (once, before parallel jobs) ---
tmp_db="/tmp/$(basename "$k2_db")"
mkdir -p "$tmp_db"
cp -r "$k2_db"/*.k2d "$tmp_db/"

mkdir -p "$outdir"

# --- Run in parallel ---
xargs -P "$parallel_jobs" -I {} bash -c 'run_kraken2 "{}"' < "$list_fq"
