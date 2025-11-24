#!/bin/bash -l
set -euo pipefail
set -x

# From launcher
SAMPLE_LIST="${SAMPLE_LIST}"       # one SAMPLE per line
FASTP_DIR="${FASTP_DIR}"           # OUT_ROOT/01_fastp
OUT_DIR="${OUT_DIR}"               # OUT_ROOT/02_prinseq
PRINSEQ_LITE="${PRINSEQ_LITE}"

PRINSEQ_COMPLEXITY_METHOD="${PRINSEQ_COMPLEXITY_METHOD:-dust}"
PRINSEQ_COMPLEXITY_THRESHOLD="${PRINSEQ_COMPLEXITY_THRESHOLD:-4}"
PRINSEQ_MIN_LEN="${PRINSEQ_MIN_LEN:-30}"
PRINSEQ_DEREP="${PRINSEQ_DEREP:-14}"

# Optional TSV mapping: SAMPLE <tab> ABS_MERGED_FASTQ
PRINSEQ_INPUTS_TSV="${PRINSEQ_INPUTS_TSV:-}"

SAMPLE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLE_LIST}")"
[[ -n "${SAMPLE}" ]] || { echo "[ERR] empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 2; }

mkdir -p "${OUT_DIR}"

# Resolve merged input fastq(.gz)
if [[ -n "${PRINSEQ_INPUTS_TSV}" && -f "${PRINSEQ_INPUTS_TSV}" ]]; then
  IN_MERGED="$(awk -F'\t' -v s="${SAMPLE}" '$1==s {print $2}' "${PRINSEQ_INPUTS_TSV}" | head -n1 || true)"
else
  IN_MERGED="${FASTP_DIR}/${SAMPLE}/${SAMPLE}_merged.fastq.gz"
fi
[[ -f "${IN_MERGED}" ]] || { echo "[ERR] merged fastq missing: ${IN_MERGED}"; exit 3; }

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

# 2) Header fix for DeDup
#gzip -df "${GOOD}.fastq.gz"
#sed -i 's/^@/@M_/;n;n;n' "${GOOD}.fastq"
#awk 'NR%4==1{sub(/^@/,"@M_")}1' "${GOOD}.fastq" > "${GOOD}.fixed.fastq"
#mv -f "${GOOD}.fixed.fastq" "${GOOD}.fastq"

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