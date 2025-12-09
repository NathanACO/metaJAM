#!/bin/bash -l
set -euo pipefail
set -x

# From launcher
FASTP_DIR="${FASTP_DIR}"           # OUT_ROOT/01_fastp
OUT_DIR="${OUT_DIR}"               # OUT_ROOT/02_prinseq
PRINSEQ_LITE="${PRINSEQ_LITE}"

PRINSEQ_COMPLEXITY_METHOD="${PRINSEQ_COMPLEXITY_METHOD:-dust}"
PRINSEQ_COMPLEXITY_THRESHOLD="${PRINSEQ_COMPLEXITY_THRESHOLD:-4}"
PRINSEQ_MIN_LEN="${PRINSEQ_MIN_LEN:-30}"
PRINSEQ_DEREP="${PRINSEQ_DEREP:-14}"

# Optional TSV mapping: SAMPLE <tab> ABS_MERGED_FASTQ
PRINSEQ_INPUTS_TSV="${PRINSEQ_INPUTS_TSV:?}"

idx="${SLURM_ARRAY_TASK_ID:-0}"
line=$(sed -n "$((idx+1))p" "${PRINSEQ_INPUTS_TSV}" || true)
if [[ -z "$line" ]]; then
  echo "[WARN] No line for array index ${idx} in ${PRINSEQ_INPUTS_TSV}, exiting."
  exit 0
fi

SAMPLE=$(printf '%s\n' "$line" | cut -f1)
IN_MERGED=$(printf '%s\n' "$line" | cut -f2-)

if [[ -z "$IN_MERGED" ]]; then
  echo "[ERR] merged path empty for sample ${SAMPLE} (line: ${line})"
  exit 3
fi

if [[ ! -f "$IN_MERGED" ]]; then
  echo "[ERR] merged fastq missing: ${IN_MERGED}"
  ls -l "$(dirname "$IN_MERGED")" || true
  exit 3
fi

# File stems (NO per-sample folder)
GOOD="${OUT_DIR}/${SAMPLE}_merged.complexity_filtered"
BAD="${OUT_DIR}/${SAMPLE}_merged.low_complexity"
DEREP_OUT="${OUT_DIR}/${SAMPLE}_merged.complexity_filtered.duplicatesremoved"

# 1) Low complexity filter
zcat "${IN_MERGED}" | perl "${PRINSEQ_LITE}" \
  -fastq stdin \
  -out_good "${GOOD}" \
  -out_bad "${BAD}" \
  -lc_method "${PRINSEQ_COMPLEXITY_METHOD}" \
  -lc_threshold "${PRINSEQ_COMPLEXITY_THRESHOLD}" \
  -min_len "${PRINSEQ_MIN_LEN}" \
  -line_width 0

gzip -f "${BAD}.fastq"

# 3) Remove exact fwd/rev duplicates
perl "${PRINSEQ_LITE}" \
  -fastq "${GOOD}.fastq" \
  -out_good "${DEREP_OUT}" \
  -out_bad null \
  -min_len "${PRINSEQ_MIN_LEN}" \
  -derep "${PRINSEQ_DEREP}" \
  -line_width 0

gzip -f "${DEREP_OUT}.fastq"

echo "[OK] PRINSEQ done for ${SAMPLE}"
