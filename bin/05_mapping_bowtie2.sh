#!/bin/bash -l
set -euo pipefail
set -x

# -------------------- env from launcher --------------------
SAMPLE_LIST="${SAMPLE_LIST}"     # list with either SAMPLE IDs or absolute FASTQ paths
INPUT_DIR="${INPUT_DIR}"         # used when list provides SAMPLE IDs (normal pipeline)
OUTPUT_DIR="${OUTPUT_DIR}"       # expected: <OUT_ROOT>/04_mapping
TMPDIR="${TMPDIR}"

# Bowtie2 index prefixes (optional; only needed for defaults / special-case checks)
PHYLONORWAY="${PHYLONORWAY:-}"
PLASTID="${PLASTID:-}"
MITO="${MITO:-}"
MAM_BIRD_FISH="${MAM_BIRD_FISH:-}"

# PhyloNorway header file (for @SQ lines); keep backwards compatibility with older variable name "HEADER"
PHYLONORWAY_HEADER="${PHYLONORWAY_HEADER:-}"

# List of DB prefixes to map (space-separated). Order matters.
MAP_DB_LIST="${MAP_DB_LIST:-}"

# Modules
ml "${BOWTIE2_MODULE}"
ml "${SAMTOOLS_MODULE}"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# -------------------- helpers --------------------
db_tag_from_prefix() {
  local p="$1"
  p="${p%/}"
  local tag
  tag="$(basename "$p")"
  # common case: bowtie2 prefixes end with "_" -> remove trailing underscores
  while [[ "$tag" == *_ ]]; do tag="${tag%_}"; done
  # if tag becomes empty (very edge case), fall back to directory name
  if [[ -z "$tag" ]]; then
    tag="$(basename "$(dirname "$p")")"
  fi
  echo "$tag"
}

# -------------------- resolve input ------------------------
LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLE_LIST}")"
[[ -n "${LINE}" ]] || { echo "[ERR] empty line for task ${SLURM_ARRAY_TASK_ID}"; exit 2; }

# If LINE is an absolute path, use it as the very first input.
# Otherwise, treat it as SAMPLE and read standard Kraken output.
if [[ "${LINE}" == /* ]]; then
  IN_UNCLAS="${LINE}"
  SAMPLE="$(basename "${IN_UNCLAS}")"
  SAMPLE="${SAMPLE%_GTDB_unclas.fastq.gz}"
  SAMPLE="${SAMPLE%_merged.dust.rmdup.fastq.gz}"
else
  SAMPLE="${LINE}"
  IN_UNCLAS="${INPUT_DIR}/${SAMPLE}_GTDB_unclas.fastq.gz"
fi

[[ -f "${IN_UNCLAS}" ]] || { echo "[ERR] Missing input FASTQ: ${IN_UNCLAS}"; exit 2; }

# -------------------- DB list default ----------------------
if [[ -z "${MAP_DB_LIST}" ]]; then
  # Default to the original 4 DBs (if set in config)
  MAP_DB_LIST="${PHYLONORWAY} ${PLASTID} ${MITO} ${MAM_BIRD_FISH}"
fi

# Convert to array (paths must not contain spaces)
read -r -a DB_PREFIXES <<< "${MAP_DB_LIST}"
if [[ "${#DB_PREFIXES[@]}" -lt 1 ]]; then
  echo "[ERR] MAP_DB_LIST is empty and no defaults available."
  exit 2
fi

LAST_DB_PREFIX="${DB_PREFIXES[$((${#DB_PREFIXES[@]}-1))]}"
LAST_DB_TAG="$(db_tag_from_prefix "${LAST_DB_PREFIX}")"

# -------------------- output layout ------------------------
# New run-specific folder:
#   04_mapping/<SAMPLE>/<SAMPLE>_<LAST_DB_TAG>/
RUN_DIR="${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_${LAST_DB_TAG}"
UNCLASS_DIR="${RUN_DIR}/unclass"
mkdir -p "${RUN_DIR}" "${UNCLASS_DIR}"

echo "[INFO] SAMPLE=${SAMPLE}"
echo "[INFO] IN_UNCLAS=${IN_UNCLAS}"
echo "[INFO] MAP_DB_LIST=${MAP_DB_LIST}"
echo "[INFO] RUN_DIR=${RUN_DIR}"

# -------------------- mapping chain ------------------------
CUR_IN="${IN_UNCLAS}"
declare -a BAM_LIST=()

for DB_PREFIX in "${DB_PREFIXES[@]}"; do
  [[ -n "${DB_PREFIX}" ]] || continue

  DB_TAG="$(db_tag_from_prefix "${DB_PREFIX}")"
  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB_TAG} alignment started for ${SAMPLE}"

  SAM_NOHDR="${RUN_DIR}/${SAMPLE}.${DB_TAG}.noheader.sam"
  SAM_ALN="${RUN_DIR}/${SAMPLE}.${DB_TAG}.aln.sam"
  SAM_OUT="${RUN_DIR}/${SAMPLE}.b2.k1000.${DB_TAG}.sam"
  BAM_OUT="${RUN_DIR}/${SAMPLE}.b2.k1000.${DB_TAG}.bam"
  UN_OUT="${UNCLASS_DIR}/${SAMPLE}_${DB_TAG}_unclass.fq"
  UN_OUT_GZ="${UN_OUT}.gz"

  bowtie2 --threads "${THREADS}"     -k 1000 -x "${DB_PREFIX}"     -U "${CUR_IN}"     --un "${UN_OUT}"     --end-to-end     -S "${SAM_NOHDR}"

  if [[ $? -ne 0 ]]; then
    echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: Bowtie2 ${DB_TAG} failed for ${SAMPLE}" >&2
    exit 1
  fi

  # gzip unclassified for next iteration
  gzip -f "${UN_OUT}"
  CUR_IN="${UN_OUT_GZ}"

  # Special-case: PhyloNorway header fix (only if DB_PREFIX matches PHYLONORWAY)
  if [[ -n "${PHYLONORWAY}" && "${DB_PREFIX}" == "${PHYLONORWAY}" ]]; then
    if [[ -z "${PHYLONORWAY_HEADER}" || ! -f "${PHYLONORWAY_HEADER}" ]]; then
      echo "[ERR] PHYLONORWAY_HEADER missing (needed for PhyloNorway header fix): ${PHYLONORWAY_HEADER}" >&2
      exit 2
    fi

    awk '$3 ~ /^AP_/{print $3}' "${SAM_NOHDR}" | sort -u > "${TMPDIR}/${SAMPLE}_ap_patterns.txt"

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
    ' "${TMPDIR}/${SAMPLE}_ap_patterns.txt" "${PHYLONORWAY_HEADER}" > "${TMPDIR}/${SAMPLE}_header_matches.txt"

    {
      cat "${TMPDIR}/${SAMPLE}_header_matches.txt"
      awk '!/^@/' "${SAM_NOHDR}"
    } > "${SAM_ALN}"

    rm -f "${TMPDIR}/${SAMPLE}_ap_patterns.txt" "${TMPDIR}/${SAMPLE}_header_matches.txt" "${SAM_NOHDR}"
    mv -f "${SAM_ALN}" "${SAM_OUT}"
    echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Header fix completed for ${SAMPLE} (${DB_TAG})"
  else
    mv -f "${SAM_NOHDR}" "${SAM_OUT}"
  fi

  # BAM (mapped reads only)
  samtools view -@ "${THREADS}" -bS -F 4 "${SAM_OUT}" > "${BAM_OUT}"
  if [[ $? -ne 0 ]]; then
    echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] ERROR: BAM conversion failed for ${SAMPLE} (${DB_TAG})" >&2
    exit 1
  fi

  rm -f "${SAM_OUT}"
  BAM_LIST+=("${BAM_OUT}")

  echo "[`TZ=Europe/Stockholm date '+%a %d %b %Y %T %Z'`] Bowtie2 ${DB_TAG} alignment completed for ${SAMPLE}"
done

# -------------------- merge & name-sort --------------------
MERGED="${RUN_DIR}/${SAMPLE}.b2.k1000.all.merged.bam"
SORTED="${RUN_DIR}/${SAMPLE}.b2.k1000.all.sorted.bam"

# If only one BAM, just copy to merged name for uniformity
if [[ "${#BAM_LIST[@]}" -eq 1 ]]; then
  cp -f "${BAM_LIST[0]}" "${MERGED}"
else
  samtools merge -@ "${THREADS}" -u "${MERGED}" "${BAM_LIST[@]}"
fi

samtools sort -@ "${THREADS}" -n   -T "${TMPDIR}/tmp_${SAMPLE}"   -o "${SORTED}"   "${MERGED}"

rm -f "${MERGED}"

echo "[OK] Mapping done for ${SAMPLE}"
echo "[OK] Final BAM: ${SORTED}"
