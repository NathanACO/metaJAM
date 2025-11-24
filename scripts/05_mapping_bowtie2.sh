#!/bin/bash -l
#SBATCH -J bowtie2_arr
#SBATCH -e /dev/null
#SBATCH -o /dev/null

set -euo pipefail
set -x

# -------------------- env from launcher --------------------
SAMPLE_LIST="${SAMPLE_LIST}"     # list with either SAMPLE IDs or absolute FASTQ paths
INPUT_DIR="${INPUT_DIR}"         # used when list provides SAMPLE IDs (normal pipeline)
OUTPUT_DIR="${OUTPUT_DIR}"
TMPDIR="${TMPDIR}"

PHYLONORWAY="${PHYLONORWAY}"
HEADER="${HEADER}"
PLASTID="${PLASTID}"
MITO="${MITO}"
MAM_BIRD_FISH="${MAM_BIRD_FISH}"

# Modules
ml "${BOWTIE2_MODULE}"
ml "${SAMTOOLS_MODULE}"

# Threads (fallback to 16)
THREADS="${SLURM_CPUS_PER_TASK:-16}"

# Ensure output layout
mkdir -p "${OUTPUT_DIR}/unclass" "${OUTPUT_DIR}/${SAMPLE}"

# -------------------- resolve input ------------------------
LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLE_LIST}")"
[[ -n "${LINE}" ]] || { echo "[ERR] empty line for task ${SLURM_ARRAY_TASK_ID}"; exit 2; }

# If LINE is an absolute path, use it as the very first input.
# Otherwise, treat it as SAMPLE and read standard Kraken output.
if [[ "${LINE}" == /* ]]; then
  IN_UNCLAS="${LINE}"
  # Derive SAMPLE for naming
  SAMPLE="$(basename "${IN_UNCLAS}")"
  SAMPLE="${SAMPLE%_GTDB_unclas.fastq.gz}"
  SAMPLE="${SAMPLE%_merged.dust.rmdup.fastq.gz}"
else
  SAMPLE="${LINE}"
  IN_UNCLAS="${INPUT_DIR}/${SAMPLE}_GTDB_unclas.fastq.gz"
fi

[[ -f "${IN_UNCLAS}" ]] || { echo "[ERR] Missing input FASTQ: ${IN_UNCLAS}"; exit 2; }

DB1="PhyloNorway"
DB2="Plastid"
DB3="Mito"
DB4="MBF"

# --------------------------- PHYLONORWAY ---------------------------
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB1} alignment started for ${SAMPLE}"

bowtie2 --threads "${THREADS}" \
  -k 1000 -x "${PHYLONORWAY}" \
  -U "${IN_UNCLAS}" \
  --end-to-end \
  -S "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_PhyloNorway_noheader.sam" \
  --un "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB1}_unclass.fq"

if [[ $? -eq 0 ]]; then
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB1} alignment completed for ${SAMPLE}"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: Bowtie2 ${DB1} alignment failed for ${SAMPLE}" >&2
  exit 1
fi

# gzip unclassified
gzip -f "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB1}_unclass.fq"
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] gzipped unclassified ${DB1} FASTQ for ${SAMPLE}"

# Header fix for PhyloNorway (add only SNs that appear in the SAM)
awk '$3 ~ /^AP_/{print $3}' "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_PhyloNorway_noheader.sam" | sort -u > "${TMPDIR}/${SAMPLE}_ap_patterns.txt"

awk '
  NR==FNR { pat[$1]=1; next }
  /^@SQ/ {
    for (i=1; i<=NF; i++) {
      if ($i ~ /^SN:/) {
        split($i,a,":");
        if (a[2] in pat && !seen[a[2]]++) print;
      }
    }
  }
' "${TMPDIR}/${SAMPLE}_ap_patterns.txt" "${HEADER}" > "${TMPDIR}/${SAMPLE}_header_matches.txt"

{
  cat "${TMPDIR}/${SAMPLE}_header_matches.txt"
  awk '!/^@/' "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_PhyloNorway_noheader.sam"
} > "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_PhyloNorway_aln.sam"

rm -f "${TMPDIR}/${SAMPLE}_ap_patterns.txt" "${TMPDIR}/${SAMPLE}_header_matches.txt" "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_PhyloNorway_noheader.sam"
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Header fix completed for ${SAMPLE}"

# BAM (mapped reads only)
samtools view -@ "${THREADS}" -bS -F 4 \
  "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_PhyloNorway_aln.sam" \
  > "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.${DB1}.bam"
if [[ $? -eq 0 ]]; then
  rm -f "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_PhyloNorway_aln.sam"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: BAM conversion failed for ${SAMPLE} (${DB1})" >&2
  exit 1
fi

# ----------------------------- PLASTID -----------------------------
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB2} alignment started for ${SAMPLE}"

bowtie2 --threads "${THREADS}" \
  -k 1000 -x "${PLASTID}" \
  -U "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB1}_unclass.fq.gz" \
  --un "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB2}_unclass.fq" \
  --end-to-end \
  -S "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB2}.sam"

if [[ $? -eq 0 ]]; then
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB2} alignment completed for ${SAMPLE}"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: Bowtie2 ${DB2} alignment failed for ${SAMPLE}" >&2
  exit 1
fi

gzip -f "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB2}_unclass.fq"
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] gzipped unclassified ${DB2} FASTQ for ${SAMPLE}"

samtools view -@ "${THREADS}" -bS -F 4 \
  "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB2}.sam" \
  > "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.${DB2}.bam"
if [[ $? -eq 0 ]]; then
  rm -f "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB2}.sam"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: BAM conversion failed for ${SAMPLE} (${DB2})" >&2
  exit 1
fi

# ------------------------------- MITO ------------------------------
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB3} alignment started for ${SAMPLE}"

bowtie2 --threads "${THREADS}" \
  -k 1000 -x "${MITO}" \
  -U "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB2}_unclass.fq.gz" \
  --un "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB3}_unclass.fq" \
  --end-to-end \
  -S "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB3}.sam"

if [[ $? -eq 0 ]]; then
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB3} alignment completed for ${SAMPLE}"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: Bowtie2 ${DB3} alignment failed for ${SAMPLE}" >&2
  exit 1
fi

gzip -f "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB3}_unclass.fq"
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] gzipped unclassified ${DB3} FASTQ for ${SAMPLE}"

samtools view -@ "${THREADS}" -bS -F 4 \
  "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB3}.sam" \
  > "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.${DB3}.bam"
if [[ $? -eq 0 ]]; then
  rm -f "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB3}.sam"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: BAM conversion failed for ${SAMPLE} (${DB3})" >&2
  exit 1
fi

# ------------------------- MAM-BIRD-FISH ---------------------------
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB4} alignment started for ${SAMPLE}"

bowtie2 --threads "${THREADS}" \
  -k 1000 -x "${MAM_BIRD_FISH}" \
  -U "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB3}_unclass.fq.gz" \
  --un "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB4}_unclass.fq" \
  --end-to-end \
  -S "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB4}.sam"

if [[ $? -eq 0 ]]; then
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB4} alignment completed for ${SAMPLE}"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: Bowtie2 ${DB4} alignment failed for ${SAMPLE}" >&2
  exit 1
fi

gzip -f "${OUTPUT_DIR}/unclass/${SAMPLE}_${DB4}_unclass.fq"
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] gzipped unclassified ${DB4} FASTQ for ${SAMPLE}"

samtools view -@ "${THREADS}" -bS -F 4 \
  "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB4}.sam" \
  > "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.${DB4}.bam"
if [[ $? -eq 0 ]]; then
  rm -f "${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.b2.k1000.${DB4}.sam"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: BAM conversion failed for ${SAMPLE} (${DB4})" >&2
  exit 1
fi

# -------------------------- MERGE + SORT ---------------------------
echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] BAM merge and sorting started for ${SAMPLE}"

# Uncompressed merge for speed
samtools merge -@ "${THREADS}" -u \
  "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.all.merged.bam" \
  "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.${DB1}.bam" \
  "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.${DB2}.bam" \
  "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.${DB3}.bam" \
  "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.${DB4}.bam"

if [[ $? -ne 0 ]]; then
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: BAM merge failed for ${SAMPLE}" >&2
  exit 1
fi

# Name-sort the merged BAM
samtools sort -@ "${THREADS}" -n \
  -T "${TMPDIR}/tmp_${SAMPLE}" \
  -o "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.all.sorted.bam" \
  "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.all.merged.bam"

if [[ $? -eq 0 ]]; then
  rm -f "${OUTPUT_DIR}/${SAMPLE}.b2.k1000.all.merged.bam"
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] BAM merge and sorting completed for ${SAMPLE}"
else
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: BAM sort failed for ${SAMPLE}" >&2
  exit 1
fi

echo "[OK] Mapping done for ${SAMPLE}"