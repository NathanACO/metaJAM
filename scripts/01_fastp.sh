#!/bin/bash -l
set -euo pipefail
set -x

ml "${FASTP_MODULE}"

sample="${SAMPLE}"
r1="${R1_IN}"
r2="${R2_IN}"
output="${OUTPUT}"

mkdir -p "${output}"

[[ -f "$r1" ]] || { echo "[ERROR] R1 not found: $r1"; exit 1; }
[[ -f "$r2" ]] || { echo "[ERROR] R2 not found: $r2"; exit 1; }

fastp -i "${r1}" -I "${r2}" -p -c --merge \
      --overlap_len_require "${FASTP_OVERLAP_LEN_REQUIRE}" \
      --overlap_diff_limit 1 \
      --merged_out="${output}/${sample}_merged.fastq.gz" \
      --include_unmerged \
      --trim_poly_x \
      -h "${output}/${sample}.html" -R "${sample} QC Report" \
      -w 10 \
      -l "${FASTP_MIN_LENGTH}"

mkdir -p "${output}/html"
mv "${output%/}"/*.html "${output%/}/html" || true
set +x
