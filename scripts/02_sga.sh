#!/bin/bash -l
set -euo pipefail
set -x

${CONDA_INIT}
conda activate "${CONDA_ENV_SGA}"

sample="${SAMPLE}"
input="${INPUT_MERGED}"    # path to <sample>_merged.fastq.gz
output="${OUTPUT_DIR}"     # output directory

mkdir -p "${output}"


# Preprocess + index + filter
sga preprocess "${input}" \
  -o "${output}/${sample}_merged.dust.fastq.gz" \
  -m 30 \
  --dust \
  --dust-threshold="${SGA_DUST_THRESHOLD}" \
  --no-primer-check \
  --min-length=30 | tee -a "${output}/sga_${sample}.log"

sga index "${output}/${sample}_merged.dust.fastq.gz" \
  --algorithm=ropebwt \
  --threads=20 \
  -p "${output}/${sample}_merged.dust"

sga filter --homopolymer-check \
  --no-kmer-check \
  --threads 20 \
  -p "${output}/${sample}_merged.dust" \
  "${output}/${sample}_merged.dust.fastq.gz" \
  -o "${output}/${sample}_merged.dust.rmdup.fastq.gz"

set +x
conda deactivate
