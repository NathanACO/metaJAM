#!/bin/bash -l

set -euo pipefail
set -x

# ---------------- env from launcher ----------------
SAMPLE_LIST="${SAMPLE_LIST}"        # one line per item: either SAMPLE or /abs/path/to/<S>_merged.dust.rmdup.*q.gz
INPUT_DIR="${INPUT_DIR}"          # <-- was SGA_DIR
INPUT_MODE="${INPUT_MODE}"
OUTDIR="${OUTPUT_DIR}"              # OUT_ROOT/03_kraken_gtdb
KRAKEN2_MODULE="${KRAKEN2_MODULE}"
PDC_MODULE="${PDC_MODULE}"
GTDB_DB="${GTDB_SRC}"               # GTDB kraken2 db root
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"  # robust default
LABEL1="GTDB"

mkdir -p "${OUTDIR}"

# Load modules (silently ok if already loaded)
ml "${PDC_MODULE}"
ml "${KRAKEN2_MODULE}"

# ---------------- main loop ----------------
while read -r LINE; do
  [[ -z "${LINE}" ]] && continue

  SAMPLE=""
  IN_FASTQ=""

  if [[ "${LINE}" == /* ]]; then
    # Absolute path to merged input
    IN_FASTQ="${LINE}"
    b="${LINE##*/}"
    SAMPLE="${b%_merged.dust.rmdup.fastq.gz}"
    SAMPLE="${SAMPLE%_merged.dust.rmdup.fq.gz}"
    SAMPLE="${SAMPLE%_merged.complexity_filtered.duplicatesremoved.fastq.gz}"
    SAMPLE="${SAMPLE%_merged.complexity_filtered.duplicatesremoved.fq.gz}"
  else
    SAMPLE="${LINE}"
    case "${INPUT_MODE}" in
      prinseq)
        IN_FASTQ="$(compgen -G "${INPUT_DIR}/${SAMPLE}_merged.complexity_filtered.duplicatesremoved.f*q.gz" | head -n1 || true)"
        ;;
      sga)
        IN_FASTQ="$(compgen -G "${INPUT_DIR}/${SAMPLE}/${SAMPLE}_merged.dust.rmdup.f*q.gz" | head -n1 || true)"
        ;;
      auto|*)
        # Try PRINSEQ first, then SGA
        IN_FASTQ="$(compgen -G "${INPUT_DIR}/${SAMPLE}_merged.complexity_filtered.duplicatesremoved.f*q.gz" | head -n1 || true)"
        if [[ -z "${IN_FASTQ}" || ! -f "${IN_FASTQ}" ]]; then
          IN_FASTQ="$(compgen -G "${INPUT_DIR}/${SAMPLE}/${SAMPLE}_merged.dust.rmdup.f*q.gz" | head -n1 || true)"
        fi
        ;;
    esac
  fi

  if [[ -z "${IN_FASTQ}" || ! -f "${IN_FASTQ}" ]]; then
    echo "[WARN] Missing input for ${SAMPLE}; skipping."
    continue
  fi

  time_start=$(date +%s)

  kraken2 --db "${GTDB_DB}" \
          --report "${OUTDIR}/${SAMPLE}_${LABEL1}_report.txt" \
          --report-minimizer-data \
          --gzip-compressed \
          --threads "${THREADS}" \
          --output "${OUTDIR}/${SAMPLE}_${LABEL1}_output.txt" \
          --classified-out "${OUTDIR}/${SAMPLE}_${LABEL1}_clas.fastq" \
          --unclassified-out "${OUTDIR}/${SAMPLE}_${LABEL1}_unclas.fastq" \
          --memory-mapping "${IN_FASTQ}"

  pigz -p "${THREADS}" -f "${OUTDIR}/${SAMPLE}"*.fastq || true
  pigz -p "${THREADS}" -f "${OUTDIR}/${SAMPLE}"*.output.txt || true

  time_end=$(date +%s)
  echo "Time spent for Kraken2 on ${SAMPLE}: $((time_end - time_start)) seconds"
done < "${SAMPLE_LIST}"

set +x
