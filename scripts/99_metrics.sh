#!/bin/bash -l
set -euo pipefail

# -------- expected env from launcher --------
OUT_ROOT="${OUT_ROOT:?}"               # pipeline root (same as used elsewhere)
LOG_ROOT="${LOG_ROOT:-${OUT_ROOT}/00_logs}"
PRIMARY_LIST_PATH="${PRIMARY_LIST_PATH:-${OUT_ROOT}/00_logs/samples.primary.txt}"

# Optional: override list used ONLY by 99_metrics.sh
# Format per line:
#   SAMPLE<TAB>/abs/path/to/SAMPLE.small.bam
#   or: /abs/path/to/SAMPLE.small.bam
OVERRIDE_LIST_METRICS="${OVERRIDE_LIST_METRICS:-}"

# Where to write metrics
OUT_DIR="${OUT_ROOT}/99_metrics"
mkdir -p "${OUT_DIR}"
TSV="${OUT_DIR}/metrics.tsv"
NEW_ROWS="${OUT_DIR}/metrics.this_run.tsv"

# -------- Helpers --------
count_fastq () {
  local f="$1"
  [[ -s "$f" ]] || { echo "NA"; return; }
  if gzip -t "$f" >/dev/null 2>&1; then
    gzip -cd "$f" 2>/dev/null | awk 'END{print (NR?NR/4:"NA")}'
  else
    awk 'END{print (NR?NR/4:"NA")}' "$f"
  fi
}

# Load samtools module if provided
if [[ -n "${SAMTOOLS_MODULE:-}" ]]; then
  ml "${SAMTOOLS_MODULE}"
fi

have_samtools () { command -v samtools >/dev/null 2>&1; }

count_bam () {
  local b="$1"
  [[ -s "$b" ]] || { echo "NA"; return; }
  if have_samtools; then
    # -F 260 = exclude unmapped (4) + secondary (256) → primary mapped alignments only
    samtools view -c -F 260 "$b" 2>/dev/null || echo "NA"
  else
    echo "NA"
  fi
}

# Grab first existing file from a list of candidates
first_file () {
  local f
  for f in "$@"; do
    [[ -s "$f" ]] && { echo "$f"; return; }
  done
  echo ""
}

# Normalize a sample name from a read path
sample_from_read_path () {
  local p="$1" b
  b="${p##*/}"
  b="${b%_R1_001.fastq.gz}"; b="${b%_R2_001.fastq.gz}"
  b="${b%_R1.fastq.gz}";     b="${b%_R2.fastq.gz}"
  b="${b%_1.fastq.gz}";      b="${b%_2.fastq.gz}"
  b="${b%_R1_001.fq.gz}";    b="${b%_R2_001.fq.gz}"
  b="${b%_R1.fq.gz}";        b="${b%_R2.fq.gz}"
  b="${b%_1.fq.gz}";         b="${b%_2.fq.gz}"
  echo "$b"
}

# Build allowed sample set from PRIMARY_LIST_PATH (List_fastq.txt)
declare -A ALLOWED=()
if [[ -s "${PRIMARY_LIST_PATH}" ]]; then
  while IFS= read -r p; do
    [[ -z "$p" ]] && continue
    s="$(sample_from_read_path "$p")"
    [[ -n "$s" ]] && ALLOWED["$s"]=1
  done < <(awk 'NF>0 && $0 !~ /^#/' "${PRIMARY_LIST_PATH}" 2>/dev/null)
fi

# -------- Build sample set (union of sources OR override) --------
declare -A SAMPLES=()
declare -A OVERRIDE_SAMPLES=()
declare -A BAMDAM_BAM_OVERRIDE=()

# If OVERRIDE_LIST_METRICS is set, we ONLY use that list
# Lines can be:
#   SAMPLE<TAB>/abs/path/to/SAMPLE.small.bam
#   /abs/path/to/SAMPLE.small.bam
if [[ -n "${OVERRIDE_LIST_METRICS}" && -f "${OVERRIDE_LIST_METRICS}" ]]; then
  while IFS=$'\t' read -r col1 col2 || [[ -n "${col1}" ]]; do
    [[ -z "${col1}" ]] && continue

    local sample bam
    if [[ -n "${col2}" ]]; then
      sample="${col1}"
      bam="${col2}"
    else
      bam="${col1}"
      b="${bam##*/}"
      sample="${b%%.*}"   # up to first dot = sample prefix
    fi

    if [[ ! -f "${bam}" ]]; then
      echo "[WARN] metrics override BAM not found: ${bam}"
      continue
    fi

    OVERRIDE_SAMPLES["${sample}"]=1
    SAMPLES["${sample}"]=1
    BAMDAM_BAM_OVERRIDE["${sample}"]="${bam}"
  done < <(awk 'NF>0 && $0 !~ /^#/' "${OVERRIDE_LIST_METRICS}")

else
  # -------- Normal mode: union of sources, but restricted to ALLOWED samples --------

  # Helper: add sample only if in ALLOWED (or if ALLOWED is empty)
  add_sample_if_allowed () {
    local s="$1"
    if [[ ${#ALLOWED[@]} -eq 0 || -n "${ALLOWED[$s]:-}" ]]; then
      SAMPLES["$s"]=1
    fi
  }

  # 1) From PRIMARY list (if present)
  if [[ -s "${PRIMARY_LIST_PATH}" ]]; then
    for s in "${!ALLOWED[@]}"; do
      SAMPLES["$s"]=1
    done
  fi

  # 2) From 01_fastp/<sample>/_merged.fastq.gz
  for f in "${OUT_ROOT}/01_fastp"/*/*_merged.f*q.gz; do
    [[ -e "$f" ]] || continue
    b="${f##*/}"; s="${b%_merged.fastq.gz}"; s="${s%_merged.fq.gz}"
    add_sample_if_allowed "$s"
  done

  # 3) From 02_sga/<sample>/_merged.dust*.fastq.gz
  for f in "${OUT_ROOT}/02_sga"/*/*_merged.dust*.f*q.gz; do
    [[ -e "$f" ]] || continue
    b="${f##*/}"; s="${b%_merged.dust.rmdup.fastq.gz}"; s="${s%_merged.dust.rmdup.fq.gz}"
    s="${s%_merged.dust.fastq.gz}"; s="${s%_merged.dust.fq.gz}"
    add_sample_if_allowed "$s"
  done

  # 4) From 02_prinseq/*_merged.complexity_filtered*.fastq.gz
  for f in "${OUT_ROOT}/02_prinseq"/*_merged.complexity_filtered*.f*q.gz; do
    [[ -e "$f" ]] || continue
    b="${f##*/}"
    s="${b%_merged.complexity_filtered.duplicatesremoved.fastq.gz}"
    s="${s%_merged.complexity_filtered.duplicatesremoved.fq.gz}"
    s="${s%_merged.complexity_filtered.fastq.gz}"
    s="${s%_merged.complexity_filtered.fq.gz}"
    add_sample_if_allowed "$s"
  done

  # 5) From 03_kraken_gtdb/<sample>_GTDB_(clas|unclas).fastq.gz
  for f in "${OUT_ROOT}/03_kraken_gtdb"/*_GTDB_*clas.f*q.gz; do
    [[ -e "$f" ]] || continue
    b="${f##*/}"; s="${b%_GTDB_clas.fastq.gz}"; s="${s%_GTDB_clas.fq.gz}"
    s="${s%_GTDB_unclas.fastq.gz}"; s="${s%_GTDB_unclas.fq.gz}"
    add_sample_if_allowed "$s"
  done

  # 6) From 04_mapping/<sample>/*.bam
  for f in "${OUT_ROOT}/04_mapping"/*/*.bam; do
    [[ -e "$f" ]] || continue
    s="$(basename "$(dirname "$f")")"
    add_sample_if_allowed "$s"
  done

  # 7) From 05_filtering/{filterBAM,bamdam}/<sample>*.bam
  for f in "${OUT_ROOT}/05_filtering/filterBAM"/*.bam "${OUT_ROOT}/05_filtering/bamdam"/*.bam; do
    [[ -e "$f" ]] || continue
    b="${f##*/}"; s="${b%%.*}"
    add_sample_if_allowed "$s"
  done
fi

# -------- Collect metrics per sample (this run) --------
# TSV header (column order is important for merging)
HEADER="sample	raw_reads	after_fastp	after_sga_preprocess	after_sga_filter	after_prinseq	kraken_classified	kraken_unclassified	mapped_phylonorway	mapped_mito	mapped_plastid	mapped_mbf	after_filterBAM	after_bamdam"

# Start a fresh per-run rows file (no header here)
: > "${NEW_ROWS}"

# Index primary list into an associative map sample->R1 path (best-effort)
declare -A R1_PATH
if [[ -s "${PRIMARY_LIST_PATH}" ]]; then
  while IFS= read -r p; do
    [[ -z "$p" ]] && continue
    s="$(sample_from_read_path "$p")"
    # Keep the first seen R1 for the sample
    if [[ -n "$s" && -z "${R1_PATH[$s]:-}" ]]; then
      R1_PATH["$s"]="$p"
    fi
  done < <(awk 'NF>0 && $0 !~ /^#/' "${PRIMARY_LIST_PATH}" 2>/dev/null)
fi

# Pick which samples to iterate:
#  - if OVERRIDE_LIST_METRICS was used -> only those
#  - else -> all discovered samples
if [[ ${#OVERRIDE_SAMPLES[@]} -gt 0 ]]; then
  sample_list=$(printf "%s\n" "${!OVERRIDE_SAMPLES[@]}" | sort)
else
  sample_list=$(printf "%s\n" "${!SAMPLES[@]}" | sort)
fi

# Walk samples in alpha order for reproducibility
for sample in $sample_list; do
  # ---- raw reads (use PRIMARY R1 if we have it) ----
  RAW_R1="${R1_PATH[$sample]:-}"
  raw_reads="NA"
  [[ -n "$RAW_R1" ]] && raw_reads="$(count_fastq "$RAW_R1")"

  # ---- after fastp ----
  f_fastp="$(first_file \
    "${OUT_ROOT}/01_fastp/${sample}/${sample}_merged.fastq.gz" \
    "${OUT_ROOT}/01_fastp/${sample}/${sample}_merged.fq.gz")"
  after_fastp="$( [[ -n "$f_fastp" ]] && count_fastq "$f_fastp" || echo NA )"

  # ---- after SGA (preprocess + filter) ----
  f_sga_pre="$(first_file \
    "${OUT_ROOT}/02_sga/${sample}/${sample}_merged.dust.fastq.gz" \
    "${OUT_ROOT}/02_sga/${sample}/${sample}_merged.dust.fq.gz")"
  after_sga_preprocess="$( [[ -n "$f_sga_pre" ]] && count_fastq "$f_sga_pre" || echo NA )"

  f_sga_filt="$(first_file \
    "${OUT_ROOT}/02_sga/${sample}/${sample}_merged.dust.rmdup.fastq.gz" \
    "${OUT_ROOT}/02_sga/${sample}/${sample}_merged.dust.rmdup.fq.gz")"
  after_sga_filter="$( [[ -n "$f_sga_filt" ]] && count_fastq "$f_sga_filt" || echo NA )"

  # ---- after PRINSEQ ----
  f_prinseq="$(first_file \
    "${OUT_ROOT}/02_prinseq/${sample}_merged.complexity_filtered.duplicatesremoved.fastq.gz" \
    "${OUT_ROOT}/02_prinseq/${sample}_merged.complexity_filtered.duplicatesremoved.fq.gz" \
    "${OUT_ROOT}/02_prinseq/${sample}_merged.complexity_filtered.fastq.gz" \
    "${OUT_ROOT}/02_prinseq/${sample}_merged.complexity_filtered.fq.gz")"
  after_prinseq="$( [[ -n "$f_prinseq" ]] && count_fastq "$f_prinseq" || echo NA )"

  # ---- Kraken GTDB ----
  f_kraken_clas="$(first_file \
    "${OUT_ROOT}/03_kraken_gtdb/${sample}_GTDB_clas.fastq.gz" \
    "${OUT_ROOT}/03_kraken_gtdb/${sample}_GTDB_clas.fq.gz")"
  f_kraken_unclas="$(first_file \
    "${OUT_ROOT}/03_kraken_gtdb/${sample}_GTDB_unclas.fastq.gz" \
    "${OUT_ROOT}/03_kraken_gtdb/${sample}_GTDB_unclas.fq.gz")"
  kraken_classified="$( [[ -n "$f_kraken_clas"   ]] && count_fastq "$f_kraken_clas"   || echo NA )"
  kraken_unclassified="$( [[ -n "$f_kraken_unclas" ]] && count_fastq "$f_kraken_unclas" || echo NA )"

  # ---- Mapping BAMs (per DB) ----
  b_pn="${OUT_ROOT}/04_mapping/${sample}/${sample}.b2.k1000.PhyloNorway.bam"
  b_mi="${OUT_ROOT}/04_mapping/${sample}/${sample}.b2.k1000.Mito.bam"
  b_pl="${OUT_ROOT}/04_mapping/${sample}/${sample}.b2.k1000.Plastid.bam"
  b_mb="${OUT_ROOT}/04_mapping/${sample}/${sample}.b2.k1000.MBF.bam"
  b_all="${OUT_ROOT}/04_mapping/${sample}/${sample}.b2.k1000.all.sorted.bam"

  mapped_phylonorway="$(count_bam "$b_pn")"
  mapped_mito="$(count_bam "$b_mi")"
  mapped_plastid="$(count_bam "$b_pl")"
  mapped_mbf="$(count_bam "$b_mb")"
  mapped_all="$(count_bam "$b_all")"

  # ---- Post-filtering ----
  b_filtered="${OUT_ROOT}/05_filtering/filterBAM/${sample}.b2.k1000.all.filtered.bam"

  # If override file provided a BAM path for this sample, use it; otherwise default
  if [[ -n "${BAMDAM_BAM_OVERRIDE[$sample]:-}" ]]; then
    b_bamdam="${BAMDAM_BAM_OVERRIDE[$sample]}"
  else
    b_bamdam="${OUT_ROOT}/05_filtering/bamdam/${sample}/${sample}.small.bam"
  fi

  after_filterBAM="$(count_bam "$b_filtered")"
  after_bamdam="$(count_bam "$b_bamdam")"

  # Emit row (NO header here)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample" \
    "$raw_reads" \
    "$after_fastp" \
    "$after_sga_preprocess" \
    "$after_sga_filter" \
    "$after_prinseq" \
    "$kraken_classified" \
    "$kraken_unclassified" \
    "$mapped_phylonorway" \
    "$mapped_mito" \
    "$mapped_plastid" \
    "$mapped_mbf" \
    "$mapped_all" \
    "$after_filterBAM" \
    "$after_bamdam" \
    >> "${NEW_ROWS}"
done

# -------- Merge with cumulative TSV (idempotent replace-per-sample) --------
HEADER_LINE="sample	raw_reads	after_fastp	after_sga_preprocess	after_sga_filter	after_prinseq	kraken_classified	kraken_unclassified	mapped_phylonorway	mapped_mito	mapped_plastid	mapped_mbf  mapped_all	after_filterBAM	after_bamdam"
if [[ ! -s "${TSV}" ]]; then
  printf '%s\n' "${HEADER_LINE}" > "${TSV}"
fi

tmp_body="$(mktemp)"
tmp_tsv="$(mktemp)"

# Combine old body (if any) with new rows; last occurrence wins
{
  tail -n +2 "${TSV}" 2>/dev/null || true
  # NEW_ROWS has no header
  cat "${NEW_ROWS}" 2>/dev/null || true
} |
awk -F'\t' '{
  row[$1]=$0
}
END{
  for (k in row) print row[k]
}' |
sort -t $'\t' -k1,1 > "${tmp_body}"

# Reattach header atomically
printf '%s\n' "${HEADER_LINE}" > "${tmp_tsv}"
cat "${tmp_body}" >> "${tmp_tsv}"
mv -f "${tmp_tsv}" "${TSV}"
rm -f "${tmp_body}"

echo "[METRICS] Updated ${TSV}"
