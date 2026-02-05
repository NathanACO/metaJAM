#!/bin/bash -l
set -euo pipefail
set -x

# -------- env from launcher --------
SAMPLE_LIST="${SAMPLE_LIST}"
MAP_DIR="${MAP_DIR}"
OUT_DIR="${OUT_DIR}"

# Which mapping run folder to use (recommended when running filtering alone).
# If empty, we'll auto-detect the most recent mapping folder for each sample.
MAP_LAST_DB_TAG="${MAP_LAST_DB_TAG:-}"

# Toggles (honor config if set; default otherwise)
ENABLE_FILTERBAM="${ENABLE_FILTERBAM:-0}"
ENABLE_NGSLCA="${ENABLE_NGSLCA:-1}"
ENABLE_BAMDAM="${ENABLE_BAMDAM:-0}"

# Conda envs
CONDA_ENV_FILTERBAM="${CONDA_ENV_FILTERBAM:-bam-filter}"
CONDA_ENV_ngsLCA="${CONDA_ENV_ngsLCA:-ngsLCA}"

# ngsLCA taxonomy
NAMES="${NAMES}"
NODES="${NODES}"
ACC2TAX="${ACC2TAX}"

# bamdam tooling + defaults (only used if ENABLE_BAMDAM=1)
THREADS="${THREADS:-16}"
BAMDAM_STRANDED="${BAMDAM_STRANDED:-ds}"     # ds/ss
BAMDAM_MINREADS="${BAMDAM_MINREADS:-5}"
BAMDAM_MAXDAMAGE="${BAMDAM_MAXDAMAGE:-0.5}"

# Used only if ENABLE_BAMDAM=1
TOP_GENUS="${TOP_GENUS:-20}"

# -------- resolve one list line: SAMPLE id OR absolute BAM path --------
LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${SAMPLE_LIST}")"
[[ -n "${LINE}" ]] || { echo "[ERR] empty line for task ${SLURM_ARRAY_TASK_ID}"; exit 2; }

# If the line is an absolute BAM path, use it as-is; otherwise treat it as SAMPLE id.
if [[ "${LINE}" == /* && "${LINE}" == *.bam ]]; then
  IN_BAM="${LINE}"
  SAMPLE="$(basename "${IN_BAM}")"
  SAMPLE="${SAMPLE%.b2.k1000.all.sorted.bam}"

  # Infer tag from folder name: .../<SAMPLE>_<TAG>/<SAMPLE>.b2.k1000.all.sorted.bam
  if [[ -z "${MAP_LAST_DB_TAG}" ]]; then
    _parent="$(basename "$(dirname "${IN_BAM}")")"
    if [[ "${_parent}" == "${SAMPLE}_"* ]]; then
      MAP_LAST_DB_TAG="${_parent#${SAMPLE}_}"
    fi
  fi
else
  SAMPLE="${LINE}"

  if [[ -n "${MAP_LAST_DB_TAG}" ]]; then
    IN_BAM="${MAP_DIR}/${SAMPLE}/${SAMPLE}_${MAP_LAST_DB_TAG}/${SAMPLE}.b2.k1000.all.sorted.bam"
  else
    # Auto-detect: choose the most recently modified merged BAM under the mapping folder
    IN_BAM="$(ls -1t "${MAP_DIR}/${SAMPLE}/${SAMPLE}_"*/"${SAMPLE}.b2.k1000.all.sorted.bam" 2>/dev/null | head -n 1 || true)"

    if [[ -n "${IN_BAM}" ]]; then
      _parent="$(basename "$(dirname "${IN_BAM}")")"
      if [[ "${_parent}" == "${SAMPLE}_"* ]]; then
        MAP_LAST_DB_TAG="${_parent#${SAMPLE}_}"
      fi
    fi
  fi
fi

[[ -n "${MAP_LAST_DB_TAG}" ]] || MAP_LAST_DB_TAG="unknown"

# Output layout (same top-level structure as before, with run-scoped subfolders)
FILTER_DIR="${OUT_DIR}/filterBAM/${SAMPLE}/${SAMPLE}_${MAP_LAST_DB_TAG}"
NGSLCA_DIR="${OUT_DIR}/ngsLCA/${SAMPLE}/${SAMPLE}_${MAP_LAST_DB_TAG}"
BAMDAM_DIR="${OUT_DIR}/bamdam/${SAMPLE}/${SAMPLE}_${MAP_LAST_DB_TAG}"

mkdir -p "${NGSLCA_DIR}"
BAM_FOR_LCA="${IN_BAM}"
LCA_PREFIX="${NGSLCA_DIR}/${SAMPLE}.b2.k1000.all.sorted"

[[ -f "${IN_BAM}" ]] || { echo "[ERR] Missing BAM: ${IN_BAM}"; exit 2; }
echo "[INFO] SAMPLE=${SAMPLE}"
echo "[INFO] TAG=${MAP_LAST_DB_TAG}"
echo "[INFO] BAM=${IN_BAM}"


# -------- filterBAM (optional) --------
if [[ "${ENABLE_FILTERBAM}" -eq 1 ]]; then
  mkdir -p "${FILTER_DIR}"
  set +u
  conda activate "${CONDA_ENV_FILTERBAM}"
  set -u
  filterBAM filter \
    --bam "${IN_BAM}" \
    --bam-filtered "${FILTER_DIR}/${SAMPLE}.b2.k1000.all.filtered.bam" \
    --stats "${FILTER_DIR}/${SAMPLE}.b2.k1000.all.stats.tsv.gz" \
    --stats-filtered "${FILTER_DIR}/${SAMPLE}.b2.k1000.all.stats-filtered.tsv.gz" \
    --threads "${THREADS}" \
    --min-read-count 3 \
    --min-read-ani 94 \
    --min-expected-breadth-ratio 0.5 \
    --min-normalized-entropy auto \
    --min-normalized-gini auto \
    --min-breadth 0 \
    --min-avg-read-ani 90 \
    --min-coverage-evenness 0.3 \
    --min-coverage-mean 0 \
    --include-low-detection \
    --sort-by-name
  BAM_FOR_LCA="${FILTER_DIR}/${SAMPLE}.b2.k1000.all.filtered.bam"
  mkdir -p "${NGSLCA_DIR}"
  LCA_PREFIX="${NGSLCA_DIR}/${SAMPLE}.b2.k1000.all.filtered"
  set +u
  conda deactivate
  set -u
fi

# -------- ngsLCA (optional) --------
if [[ "${ENABLE_NGSLCA}" -eq 1 ]]; then
  mkdir -p "${NGSLCA_DIR}"
  set +u
  conda activate "${CONDA_ENV_ngsLCA}"
  set -u
  ngsLCA \
    -simscorelow 0.95 \
    -simscorehigh 1 \
    -names "${NAMES}" \
    -nodes "${NODES}" \
    -acc2tax "${ACC2TAX}" \
    -bam "${BAM_FOR_LCA}" \
    -outnames "${LCA_PREFIX}" \
    -fix_ncbi 0
  set +u
  conda deactivate
  set -u
fi

# -------- bamdam full per-sample pipeline (shrink + compute + krona) --------
if [[ "${ENABLE_BAMDAM}" -eq 1 ]]; then
  mkdir -p "${BAMDAM_DIR}"
  # Needs a matching LCA file
  IN_LCA="${LCA_PREFIX}.lca"
  if [[ ! -f "${IN_LCA}" ]]; then
    echo "[WARN] ENABLE_BAMDAM=1 but missing LCA (${IN_LCA}); skipping bamdam for ${SAMPLE}"
  else
    # Load Python module (optional) and activate bamdam venv
    BAMDAM_PYTHON_MODULE="${BAMDAM_PYTHON_MODULE:-python}"
    BAMDAM_VENV="${BAMDAM_VENV:-}"
    if [[ -n "${BAMDAM_PYTHON_MODULE}" ]]; then
      ml "${BAMDAM_PYTHON_MODULE}"
    fi
    if [[ -n "${BAMDAM_VENV}" && -f "${BAMDAM_VENV}/bin/activate" ]]; then
      source "${BAMDAM_VENV}/bin/activate"
    else
      echo "[WARN] BAMDAM_VENV not set or invalid; expecting 'bamdam' in PATH"
    fi

    ml "${KRONATOOLS_MODULE:-kronatools/2.8.1}"

    # shrink
    bamdam shrink \
      --in_bam "${BAM_FOR_LCA}" \
      --in_lca "${IN_LCA}" \
      --out_bam "${BAMDAM_DIR}/${SAMPLE}.small.bam" \
      --out_lca "${BAMDAM_DIR}/${SAMPLE}.small.lca" \
      --stranded "${BAMDAM_STRANDED}" \
      --show_progress

    # compute
    bamdam compute \
      --in_bam "${BAMDAM_DIR}/${SAMPLE}.small.bam" \
      --in_lca "${BAMDAM_DIR}/${SAMPLE}.small.lca" \
      --out_tsv "${BAMDAM_DIR}/${SAMPLE}.tsv" \
      --out_subs "${BAMDAM_DIR}/${SAMPLE}.subs.txt" \
      --stranded "${BAMDAM_STRANDED}" \
      --show_progress

    # krona + HTML (per-sample)
    bamdam krona \
      --in_tsv "${BAMDAM_DIR}/${SAMPLE}.tsv" \
      --out_xml "${BAMDAM_DIR}/${SAMPLE}.xml" \
      --minreads "${BAMDAM_MINREADS}" \
      --maxdamage "${BAMDAM_MAXDAMAGE}"

    ktImportXML -o "${BAMDAM_DIR}/${SAMPLE}.html" "${BAMDAM_DIR}/${SAMPLE}.xml"

    while read -r count taxid genus; do
      bamdam plotdamage \
      --in_subs "${BAMDAM_DIR}/${SAMPLE}.subs.txt" \
      --tax "$taxid" \
      --outplot "${BAMDAM_DIR}/damageplot_${genus}_${taxid}_${SAMPLE}.png" \
      || echo "Warning: bamdam failed for taxid $taxid (continuing)"
    done < <(
    awk '{for(i=2;i<=NF;i++){n=split($i,a,":");if(n>=3 && a[3]=="genus"){id=a[1];counts[id]++;name[id]=a[2];break}}}
       END{for(id in counts)printf "%d\t%s\t%s\n",counts[id],id,name[id]}' \
    "${BAMDAM_DIR}/${SAMPLE}.small.lca" \
    | sort -nrk1,1 \
    | head -n "$TOP_GENUS"
    )
  fi
fi

echo "[OK] Filtering done for ${SAMPLE}"