#!/bin/bash -l
set -euo pipefail
set -x

# -----------------------------------------------------------------------------
# Resolve sample from SLURM array (Run_pipeline submits --array=0-(N-1))
# -----------------------------------------------------------------------------
: "${SAMPLE_LIST:?SAMPLE_LIST not set}"
sample="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLE_LIST}" || true)"
if [[ -z "${sample}" ]]; then
  echo "[ERROR] Could not resolve sample for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID} from ${SAMPLE_LIST}"
  exit 1
fi

# -----------------------------------------------------------------------------
# Locate bamdam LCA input directory for this sample
# Expected structure:
#   ${OUT_ROOT}/05_filtering/bamdam/<sample>/<sample>_<tag>/<sample>.small.lca
# -----------------------------------------------------------------------------
BAMDAM_SAMPLE_DIR="${OUT_ROOT}/05_filtering/bamdam/${sample}"
if [[ ! -d "${BAMDAM_SAMPLE_DIR}" ]]; then
  echo "[ERROR] Missing bamdam dir: ${BAMDAM_SAMPLE_DIR}"
  exit 1
fi

RUN_DIR="$(find "${BAMDAM_SAMPLE_DIR}" -maxdepth 1 -type d -name "${sample}_*" | sort | head -n 1 || true)"
if [[ -z "${RUN_DIR}" ]]; then
  echo "[ERROR] Could not find run dir ${BAMDAM_SAMPLE_DIR}/${sample}_*"
  exit 1
fi


MAP_LAST_DB_TAG="${MAP_LAST_DB_TAG}"
INPUT_DIR="${RUN_DIR}"
OUTPUT_DIR="${OUT_ROOT}/05_filtering/mmseq2/${sample}/${sample}_${MAP_LAST_DB_TAG}"
mkdir -p "${OUTPUT_DIR}"

# TMP (per-sample to reduce collisions)
TMPDIR="${TMP_ROOT}/mmseqs2/${sample}"
mkdir -p "${TMPDIR}"

# -----------------------------------------------------------------------------
# Modules
# -----------------------------------------------------------------------------
USER_MMSEQS2_DB="${MMSEQS2_DB}"
: "${MMSEQS2_DB:?MMSEQS2_DB must point to the MMseqs2 target DB (provided by module or exported)}"

ml bioinfo-tools
if [[ -n "${MMSEQS2_MODULE:-}" ]]; then ml "${MMSEQS2_MODULE}"; else ml MMseqs2/15-6f452; fi
if [[ -n "${MMSEQS2_DATA_MODULE:-}" ]]; then ml "${MMSEQS2_DATA_MODULE}"; else ml MMseqs2_data/latest; fi


# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
TOP_GENERA="${MMSEQS2_TOP_GENERA:-10}"
MAX_READS="${MMSEQS2_MAX_READS:-100}"
MIN_READS="${MMSEQS2_MIN_READS:-30}"
SEED="${MMSEQS2_SEED:-}"

THREADS="${MMSEQS2_THREADS:-7}"
MAX_SEQS="${MMSEQS2_MAX_SEQS:-300}"
MIN_LEN="${MMSEQS2_MIN_LENGTH:-30}"
S_PARAM="${MMSEQS2_S:-7.5}"
SKM="${MMSEQS2_SPACED_KMER_MODE:-1}"
SPLIT_MEM="${MMSEQS2_SPLIT_MEM_LIMIT:-220G}"

DATABASE="nt_query"

LCA_FILE="${INPUT_DIR}/${sample}.small.lca"
if [[ ! -f "${LCA_FILE}" ]]; then
  echo "[ERROR] Missing LCA file: ${LCA_FILE}"
  exit 1
fi

LCA_TSV="${INPUT_DIR}/${sample}.tsv"
if [[ ! -f "${LCA_TSV}" ]]; then
  echo "[ERROR] Missing LCA TSV file to perform damage filtering: ${LCA_TSV}"
  exit 1
fi

# -----------------------------------------------------------------------------
# 1) Build combined FASTA + expected map from LCA (top10 genus / kingdom)
# -----------------------------------------------------------------------------
SEED_ARGS=()
if [[ -n "${SEED}" ]]; then
  SEED_ARGS=(--seed "${SEED}")
fi

# Prefer genera-file if provided
GENUS_ARGS=()
if [[ -n "${MMSEQS2_GENERA_FILE:-}" && -f "${MMSEQS2_GENERA_FILE}" ]]; then
  GENUS_ARGS+=( --genera-file "${MMSEQS2_GENERA_FILE}" )
fi

#If a genera override is provided, force re-run of search/convertalis by removing prior MMseqs outputs
if (( ${#GENUS_ARGS[@]} > 0 )); then
  echo "[INFO] MMSEQS2_GENERA_FILE set; removing existing MMseqs resultDB outputs to force re-run."
  rm -f "${OUTPUT_DIR}/${DATABASE}.resultDB"*
fi

python3 "${SCRIPTS_DIR}/filter_lca_top10_subset_reads.py" \
  --lca "${LCA_FILE}" \
  --bamdam-tsv "${LCA_TSV}" \
  --min-dmg "${MMSEQS2_MIN_DMG}" \
  --top-genera "${TOP_GENERA}" \
  --max-reads "${MAX_READS}" \
  --min-reads "${MIN_READS}" \
  "${SEED_ARGS[@]}" \
  "${GENUS_ARGS[@]}" \
  --out-fasta "${OUTPUT_DIR}/${sample}.mmseqs_queries.fa" \
  --out-map   "${OUTPUT_DIR}/${sample}.mmseqs_expected.tsv"

# If nothing was written (all taxa skipped), stop gracefully
if [[ ! -s "${OUTPUT_DIR}/${sample}.mmseqs_queries.fa" ]]; then
  echo "[INFO] No queries generated for ${sample} (likely all taxa < min reads)."
  exit 0
fi

# -----------------------------------------------------------------------------
# 2) MMseqs2 search
# -----------------------------------------------------------------------------
mmseqs createdb "${OUTPUT_DIR}/${sample}.mmseqs_queries.fa" "${OUTPUT_DIR}/${DATABASE}"

if [[ -f "${OUTPUT_DIR}/${DATABASE}.resultDB.dbtype" ]]; then
  echo "[INFO] MMseqs search output exists, reusing: ${OUTPUT_DIR}/${DATABASE}.resultDB"
else
  mmseqs search \
    "${OUTPUT_DIR}/${DATABASE}" \
    "${USER_MMSEQS2_DB}" \
    "${OUTPUT_DIR}/${DATABASE}.resultDB" \
    "${TMPDIR}" \
    --threads "${THREADS}" \
    -a \
    --search-type 3 \
    --spaced-kmer-mode "${SKM}" \
    -s "${S_PARAM}" \
    --min-length "${MIN_LEN}" \
    --min-seq-id "${MIN_SEQID}" \
    --max-seqs "${MAX_SEQS}" \
    --split-memory-limit "${SPLIT_MEM}"
fi

if [[ -s "${OUTPUT_DIR}/${DATABASE}.resultDB.m8" ]]; then
  echo "[INFO] convertalis output exists, reusing: ${OUTPUT_DIR}/${DATABASE}.resultDB.m8"
else
  mmseqs convertalis \
    --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qlen,tcov,qcov,qaln,taln" \
    "${OUTPUT_DIR}/${DATABASE}" \
    "${USER_MMSEQS2_DB}" \
    "${OUTPUT_DIR}/${DATABASE}.resultDB" \
    "${OUTPUT_DIR}/${DATABASE}.resultDB.m8" \
    --format-mode 4 \
    --search-type 3
fi

## Bits filtering ##
awk -F'\t' 'NR==1 || $12 > "${MIN_BITS}"' "${OUTPUT_DIR}/${DATABASE}.resultDB.m8" > "${OUTPUT_DIR}/${DATABASE}.resultDB_bits.m8"

# -----------------------------------------------------------------------------
# 3) Add taxonomic info to MMseqs2 m8 (optional; requires TAXADB)
# -----------------------------------------------------------------------------
if [[ -n "${TAXADB_VENV:-}" && -n "${TAXADB_SQLITE:-}" && -f "${TAXADB_SQLITE}" ]]; then
  source "${TAXADB_VENV}/bin/activate"
  python3 "${SCRIPTS_DIR}/add_taxid_info.py" \
    -blst Mmseqs2 \
    -b "${OUTPUT_DIR}/${DATABASE}.resultDB_bits.m8" \
    -d "${TAXADB_SQLITE}" \
    -o "${OUTPUT_DIR}/${sample}.${DATABASE}.resultDB_bits_nucl_taxid.m8" \
    -del tab
  deactivate
else
  echo "[INFO] TAXADB not configured; skipping taxid augmentation."
fi

rm -f "${OUTPUT_DIR}/nt_query.resultDB."[0-9]* # Remove intermediate files

# -----------------------------------------------------------------------------
# Post-processing: pick one best hit per query + compare taxonomy vs expected
# -----------------------------------------------------------------------------
AMBIG_FRAC="${MMSEQS2_AMBIG_FRAC:-0.005}"

python3 "${SCRIPTS_DIR}/mmseqs_pick_besthit_with_ambiguity.py" \
  --m8 "${OUTPUT_DIR}/${sample}.${DATABASE}.resultDB_bits_nucl_taxid.m8" \
  --ambig-frac "${AMBIG_FRAC}" \
  --out "${OUTPUT_DIR}/${sample}.${DATABASE}.besthit.assigned.tsv"

python3 "${SCRIPTS_DIR}/postprocess_mmseqs_taxonomy_compare.py" \
  --expected-map "${OUTPUT_DIR}/${sample}.mmseqs_expected.tsv" \
  --assigned "${OUTPUT_DIR}/${sample}.${DATABASE}.besthit.assigned.tsv" \
  --out-summary "${OUTPUT_DIR}/${sample}.${DATABASE}.bamdam_mmseqs.evaluation.summary.tsv"

echo "[mmseqs2] Evaluation summary written: "${OUTPUT_DIR}/${sample}.${DATABASE}.bamdam_mmseqs.evaluation.summary.tsv""
