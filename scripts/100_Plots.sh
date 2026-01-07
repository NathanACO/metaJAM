#!/bin/bash -l
set -euo pipefail

# -------- expected env from launcher --------
OUT_ROOT="${OUT_ROOT:?}"               # pipeline root
LOG_ROOT="${LOG_ROOT:-${OUT_ROOT}/00_logs}"

# Default sample list path (can be overridden by env)
PRIMARY_LIST_PATH="${PRIMARY_LIST_PATH:-${OUT_ROOT}/00_Samples_prefix/samples.primary.txt}"

# Inputs produced by previous steps (can be overridden)
METRICS_TSV="${METRICS_TSV:-${OUT_ROOT}/99_metrics/metrics.tsv}"

# Optional: mapping/filtering run tag (used to select the correct run folders)
MAP_LAST_DB_TAG="${MAP_LAST_DB_TAG:-}"

BAMDAM_DIR="${BAMDAM_DIR:-${OUT_ROOT}/05_filtering/bamdam}"
MIN_READS="${PLOTS_BAMDAM_MIN_READS:-1}"    # Minimum reads per sample to include in bamdam plots
PLOTS_BAMDAM_PLOT_MODE="${PLOTS_BAMDAM_PLOT_MODE:-heatmap}"
PLOTS_DAMAGE_THRESHOLD="${PLOTS_DAMAGE_THRESHOLD:-5}"
# User-provided metadata file
METADATA_PATH="${METADATA_PATH:?Need METADATA_PATH pointing to metadata TSV}"

# Plot output dir
if [[ -n "${MAP_LAST_DB_TAG}" ]]; then
  PLOT_DIR="${PLOT_DIR:-${OUT_ROOT}/100_plots/${MAP_LAST_DB_TAG}}"
else
  PLOT_DIR="${PLOT_DIR:-${OUT_ROOT}/100_plots}"
fi

# Optional conda env
CONDA_ENV_PLOTS="${CONDA_ENV_PLOTS:-plots}"
RSCRIPT="${RSCRIPT:-Rscript}"

mkdir -p "${PLOT_DIR}" "${LOG_ROOT}"

if [[ -n "${SCRIPTS_DIR:-}" ]]; then
  BASE_SCRIPTS_DIR="${SCRIPTS_DIR}"
else
  BASE_SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

PLOTS_R="${BASE_SCRIPTS_DIR}/100_Plots.R"   # OR 100_plots.R — see note below

if [[ ! -f "${PLOTS_R}" ]]; then
  echo "[100_Plots] ERROR: cannot find plots R script at ${PLOTS_R}" >&2
  echo "[100_Plots] Contents of BASE_SCRIPTS_DIR:" >&2
  ls -l "${BASE_SCRIPTS_DIR}" >&2 || true
  exit 1
fi

echo "[100_Plots] using R script: ${PLOTS_R}"


# --- activate conda env if available ---
if [[ -n "${CONDA_ENV_PLOTS}" ]]; then
    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)"
        conda activate "${CONDA_ENV_PLOTS}" || \
          echo "[100_Plots] WARNING: could not activate env ${CONDA_ENV_PLOTS}, using system R" >&2
    else
        echo "[100_Plots] WARNING: conda not found, using system R" >&2
    fi
fi

# --- sanity checks (but do not hard-fail on missing bamdam etc.) ---
if [[ ! -f "${METRICS_TSV}" ]]; then
    echo "[100_Plots] WARNING: metrics file not found at ${METRICS_TSV}, skipping metrics plot" >&2
fi

if [[ ! -d "${BAMDAM_DIR}" ]]; then
    echo "[100_Plots] WARNING: bamdam directory not found at ${BAMDAM_DIR}, skipping heatmap" >&2
fi

# --- run the R plotting script ---
"${RSCRIPT}" "${PLOTS_R}" \
  --metrics "${METRICS_TSV}" \
  --samples "${PRIMARY_LIST_PATH}" \
  --bamdam_dir "${BAMDAM_DIR}" \
  --metadata "${METADATA_PATH}" \
  --db_tag "${MAP_LAST_DB_TAG}" \
  --outdir "${PLOT_DIR}" \
  --min_reads "${MIN_READS}" \
  --bamdam_plot "${PLOTS_BAMDAM_PLOT_MODE}" \
  --damage_threshold "${PLOTS_DAMAGE_THRESHOLD}" \
  --plot_low_damage_taxa "${PLOTS_PLOT_LOW_DAMAGE_TAXA:-1}" \
  --exclude_taxa "${PLOTS_EXCLUDE_TAXA:-}" \
  --taxa_per_plot "${BAMDAM_TAXA_PER_PLOT:-30}"
  > "${LOG_ROOT}/100_plots${MAP_LAST_DB_TAG:+.${MAP_LAST_DB_TAG}}.R.out" 2>&1

echo "[100_Plots] finished. Plots in ${PLOT_DIR}"
