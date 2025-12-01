#!/usr/bin/env bash
# =============================================================================
# Pipeline_sedaDNA — master launcher
# - Reads config.sh
# - Builds/uses a persistent PRIMARY_LIST_PATH with canonical sample IDs
# - Auto-discovers inputs for each step (fastp → SGA → Kraken2 → Mapping → Filter/ngsLCA)
# - Skips work if outputs exist unless FORCE_* = 1
# - Per-step SBATCH options (mem in G), passed via CLI (override headers)
# - Per-sample job names for FASTP & SGA; array for mapping; single jobs for Kraken & Filtering
# - Submits SLURM jobs with proper afterok dependencies, that means for instance for a sample S1, that SGA(S1) cannot start until FASTP(S1) finishes successfully.
# =============================================================================
set -euo pipefail

# 1) Take the config path as positional arg 1, then remove it from "$@"
CFG=${1:?Usage: bash run_pipeline.sh /full/path/to/config.sh}
shift

# 2) Load the config
source "$CFG"
declare -A JID_FASTP=()
declare -A JID_SGA=()
JOB_PRINSEQ=""
JID_KRAKEN=""
JOB_MAP=""
JOB_FILTER=""

# Resolve where the step scripts are
RUN_DIR="$(cd "$(dirname "$0")" && pwd)"
: "${SCRIPTS_DIR:=${RUN_DIR}/scripts}"   # fallback if not set in config

# tiny helpers
require_dir()  { [[ -d "$1" ]] || { echo "[ERROR] Missing dir: $1"; exit 1; }; }
require_file() { [[ -f "$1" ]] || { echo "[ERROR] Missing file: $1"; exit 1; }; }
require_dir "${SCRIPTS_DIR}"

# 3) (Optional but helpful) guard + echo
require_nonempty() { local k="$1"; [[ -n "${!k:-}" ]] || { echo "[CONFIG] $k is empty or unset (from: $CFG)"; exit 1; }; }
require_nonempty OUT_ROOT
require_nonempty TMP_ROOT
require_nonempty LOG_ROOT
echo "[CONFIG] using $CFG"
echo "[CONFIG] OUT_ROOT=$OUT_ROOT"
echo "[CONFIG] TMP_ROOT=$TMP_ROOT"
echo "[CONFIG] LOG_ROOT=$LOG_ROOT"

# Optional CLI overrides (you can ignore and just edit config.sh)
PERSIST_OVERRIDES=0
EXTRA_ARGS=("${@:2}")
i=0
while [[ $i -lt ${#EXTRA_ARGS[@]} ]]; do
  key="${EXTRA_ARGS[$i]}"
  case "$key" in
    --fastp-overlap) FASTP_OVERLAP_LEN_REQUIRE="${EXTRA_ARGS[$((i+1))]}"; i=$((i+2));;
    --fastp-minlen)  FASTP_MIN_LENGTH="${EXTRA_ARGS[$((i+1))]}";          i=$((i+2));;
    --sga-dust)      SGA_DUST_THRESHOLD="${EXTRA_ARGS[$((i+1))]}";        i=$((i+2));;
    --persist)       PERSIST_OVERRIDES=1;                                 i=$((i+1));;
    *) echo "[WARN] Unknown flag '$key' ignored";                         i=$((i+1));;
  esac
done

persist_kv() {
  local file="$1" key="$2" val="$3"
  if grep -qE "^[[:space:]]*${key}=" "$file"; then
    sed -i -E "s|^[[:space:]]*${key}=.*|${key}=${val}|" "$file"
  else
    echo "${key}=${val}" >> "$file"
  fi
}
if [[ ${PERSIST_OVERRIDES} -eq 1 ]]; then
  echo "[INFO] Persisting overrides into ${CFG}"
  persist_kv "$CFG" FASTP_OVERLAP_LEN_REQUIRE "${FASTP_OVERLAP_LEN_REQUIRE}"
  persist_kv "$CFG" FASTP_MIN_LENGTH "${FASTP_MIN_LENGTH}"
  persist_kv "$CFG" SGA_DUST_THRESHOLD "${SGA_DUST_THRESHOLD}"
fi

# -----------------------------------------------------------------------------
# Directories & logs
# -----------------------------------------------------------------------------
mkdir -p "${OUT_ROOT}"/{00_logs,01_fastp,02_sga,03_kraken_gtdb,04_mapping,05_filtering,99_metrics} \
         "${TMP_ROOT}" \
         "${LOG_ROOT}"/{out,error}
PRIMARY_LIST_PATH="${OUT_ROOT}/00_logs/samples.primary.txt"

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
join_by() { local IFS="$1"; shift; echo "$*"; }

sbatch_opts() {
  local step="$1"
  local acct="${SBATCH_DEFAULT_ACCOUNT}"
  local part="${SBATCH_DEFAULT_PARTITION}"
  local cpus="${SBATCH_DEFAULT_CPUS}"
  local mem="${SBATCH_DEFAULT_MEM}"      # e.g. "64G"
  local time="${SBATCH_DEFAULT_TIME}"
  local qos="${SBATCH_DEFAULT_QOS}"
  local extra="${SBATCH_DEFAULT_EXTRA}"
  local jname=""

  for v in ACCOUNT PARTITION CPUS MEM TIME QOS EXTRA JOB_NAME; do
    local var="${step}_SBATCH_${v}"
    if [[ -n "${!var:-}" ]]; then
      case "$v" in
        ACCOUNT) acct="${!var}";;
        PARTITION) part="${!var}";;
        CPUS) cpus="${!var}";;
        MEM) mem="${!var}";;
        TIME) time="${!var}";;
        QOS) qos="${!var}";;
        EXTRA) extra="${!var}";;
        JOB_NAME) jname="${!var}";;
      esac
    fi
  done

  SBATCH_BUILT_OPTS=()
  [[ -n "$acct" ]]  && SBATCH_BUILT_OPTS+=(--account="$acct")
  [[ -n "$part" ]]  && SBATCH_BUILT_OPTS+=(-p "$part")
  [[ -n "$cpus" ]]  && SBATCH_BUILT_OPTS+=(--cpus-per-task="$cpus")
  [[ -n "$mem"  ]]  && SBATCH_BUILT_OPTS+=(--mem="$mem")
  [[ -n "$time" ]]  && SBATCH_BUILT_OPTS+=(-t "$time")
  [[ -n "$qos"  ]]  && SBATCH_BUILT_OPTS+=(--qos="$qos")
  [[ -n "$jname" ]] && SBATCH_BUILT_OPTS+=(--job-name="$jname")
  if [[ -n "$extra" ]]; then
    EXTRAS_ARR=($extra)  # split on spaces
    SBATCH_BUILT_OPTS+=("${EXTRAS_ARR[@]}")
  fi
}

build_primary_list() {
  : > "${PRIMARY_LIST_PATH}"
  # 0) If user provided a list file of FASTQs, prefer it.
  #    - Accept absolute paths
  #    - Keep only R1/_1 files
  #    - Allow comments and blank lines
  if [[ -n "${READS_LIST:-}" && -f "${READS_LIST}" ]]; then
    awk 'NF>0 && $0 !~ /^#/' "${READS_LIST}" \
    | sed 's/[[:space:]]*$//' \
    | grep -E '(^/).*(_R1(_001)?|_1)\.(fastq|fq)\.gz$' \
    | while IFS= read -r p; do
        [[ -f "$p" ]] && printf '%s\n' "$p"
      done >> "${PRIMARY_LIST_PATH}" || true
  elif [[ -n "${READS_GLOB:-}" ]]; then
    # Write FULL PATHS for R1 files only (handles _R1, _R1_001, _1)
    compgen -G "${READS_GLOB}" 2>/dev/null | grep -E '(_R1(_001)?|_1)\.(fastq|fq)\.gz$' \
      >> "${PRIMARY_LIST_PATH}" || true
  elif [[ -n "${FASTP_MERGED_DIR:-}" ]];  then
    ls -1 "${FASTP_MERGED_DIR}"/*_merged.fastq.gz 2>/dev/null \
      | sed -E 's/_merged\.fastq\.gz$//' >> "${PRIMARY_LIST_PATH}" || true
  elif [[ -n "${SGA_OUT_DIR:-}" ]];       then
    ls -1 "${SGA_OUT_DIR}"/*_merged.dust.rmdup.fastq.gz 2>/dev/null \
      | sed -E 's/_merged\.dust\.rmdup\.fastq\.gz$//' >> "${PRIMARY_LIST_PATH}" || true
  elif [[ -n "${KRAKEN_UNCLAS_DIR:-}" ]]; then
    ls -1 "${KRAKEN_UNCLAS_DIR}"/*_GTDB_unclas.fastq.gz 2>/dev/null \
      | sed -E 's/_GTDB_unclas\.fastq\.gz$//' >> "${PRIMARY_LIST_PATH}" || true
  elif [[ -n "${MAPPING_BAM_DIR:-}" ]];   then
    ls -1 "${MAPPING_BAM_DIR}"/*.b2.k1000.all.sorted.bam 2>/dev/null \
      | sed -E 's/\.b2\.k1000\.all\.sorted\.bam$//' >> "${PRIMARY_LIST_PATH}" || true
  fi
  sort -u -o "${PRIMARY_LIST_PATH}" "${PRIMARY_LIST_PATH}" || true
  echo "[INFO] Primary samples: $(wc -l < "${PRIMARY_LIST_PATH}" | tr -d ' ') -> ${PRIMARY_LIST_PATH}"
}

# --- SGA: build inputs (SAMPLE \t merged.fastq.gz) ---
SGA_INPUTS="${OUT_ROOT}/00_logs/sga.inputs.tsv"
SGA_LIST="${OUT_ROOT}/00_logs/samples.for_sga.txt"

build_sga_inputs() {
  : > "${SGA_INPUTS}"

  if [[ -n "${OVERRIDE_LIST_SGA:-}" && -f "${OVERRIDE_LIST_SGA}" ]]; then
    while IFS=$'\t' read -r a b || [[ -n "$a" ]]; do
      [[ -z "$a" ]] && continue
      if [[ -n "$b" ]]; then
        s="$a"; f="$b"
      else
        f="$a"
        s="${f##*/}"; s="${s%_merged.fastq.gz}"; s="${s%_merged.fq.gz}"
      fi
      [[ -f "$f" ]] || { echo "[WARN] Missing merged file: $f"; continue; }
      printf "%s\t%s\n" "$s" "$f" >> "${SGA_INPUTS}"
    done < "${OVERRIDE_LIST_SGA}"

  else
    # --- normalize BOTH R1 and R2 from PRIMARY_LIST_PATH to <sample> ---
    declare -A allow=()
    while read -r r; do
      [[ -z "$r" ]] && continue
      s="${r##*/}"
      # strip any of: _R1/_R2 with/without _001, and _1/_2; both fastq/fq
      s="${s%_R1_001.fastq.gz}"; s="${s%_R2_001.fastq.gz}"
      s="${s%_R1.fastq.gz}";     s="${s%_R2.fastq.gz}"
      s="${s%_1.fastq.gz}";      s="${s%_2.fastq.gz}"
      s="${s%_R1_001.fq.gz}";    s="${s%_R2_001.fq.gz}"
      s="${s%_R1.fq.gz}";        s="${s%_R2.fq.gz}"
      s="${s%_1.fq.gz}";         s="${s%_2.fq.gz}"
      allow["$s"]=1
    done < "${PRIMARY_LIST_PATH}"

    # discover merged outputs from fastp and keep only allowed samples
    find "${OUT_ROOT}/01_fastp" -type f \
      \( -name '*_merged.fastq.gz' -o -name '*_merged.fq.gz' \) 2>/dev/null \
      | sort \
      | while read -r f; do
          s="${f##*/}"; s="${s%_merged.fastq.gz}"; s="${s%_merged.fq.gz}"
          [[ -n "${allow[$s]:-}" ]] || continue
          printf "%s\t%s\n" "$s" "$f" >> "${SGA_INPUTS}"
        done
  fi

  cut -f1 "${SGA_INPUTS}" | sort -u > "${SGA_LIST}" || :
  echo "[INFO] SGA inputs: $(wc -l < "${SGA_LIST}" | tr -d ' ') samples -> ${SGA_LIST}"
}

# ------------------------------- PRINSEQ inputs -------------------------------
PRINSEQ_INPUTS="${OUT_ROOT}/00_logs/prinseq.inputs.tsv"   # SAMPLE \t ABS_MERGED_FASTQ
PRINSEQ_LIST="${OUT_ROOT}/00_logs/samples.for_prinseq.txt"

build_prinseq_inputs() {
  : > "${PRINSEQ_INPUTS}"

  if [[ -n "${OVERRIDE_LIST_PRINSEQ:-}" && -f "${OVERRIDE_LIST_PRINSEQ}" ]]; then
    # Lines can be either:
    #   SAMPLE \t /abs/path/to/SAMPLE_merged.fastq.gz
    #   /abs/path/to/SAMPLE_merged.fastq.gz
    while IFS=$'\t' read -r col1 col2 || [[ -n "${col1}" ]]; do
      [[ -z "${col1}" ]] && continue

      if [[ -n "${col2}" ]]; then
        sample_id="${col1}"
        merged_path="${col2}"
      else
        merged_path="${col1}"
        # derive SAMPLE from basename of the path
        sample_id="${merged_path##*/}"
        sample_id="${sample_id%_merged.fastq.gz}"
        sample_id="${sample_id%_merged.fq.gz}"
      fi

      [[ -f "${merged_path}" ]] || { echo "[WARN] Missing merged: ${merged_path}"; continue; }
      printf "%s\t%s\n" "${sample_id}" "${merged_path}" >> "${PRINSEQ_INPUTS}"
    done < "${OVERRIDE_LIST_PRINSEQ}"

  else
    # --- normalize BOTH R1 and R2 from PRIMARY_LIST_PATH to <sample> ---
    declare -A allow=()
    while read -r r; do
      [[ -z "$r" ]] && continue
      s="${r##*/}"
      # strip any of: _R1/_R2 (with/without _001) and _1/_2; fastq or fq
      s="${s%_R1_001.fastq.gz}"; s="${s%_R2_001.fastq.gz}"
      s="${s%_R1.fastq.gz}";     s="${s%_R2.fastq.gz}"
      s="${s%_1.fastq.gz}";      s="${s%_2.fastq.gz}"
      s="${s%_R1_001.fq.gz}";    s="${s%_R2_001.fq.gz}"
      s="${s%_R1.fq.gz}";        s="${s%_R2.fq.gz}"
      s="${s%_1.fq.gz}";         s="${s%_2.fq.gz}"
      allow["$s"]=1
    done < "${PRIMARY_LIST_PATH}"

    # discover merged outputs from FASTP and keep only allowed samples
    find "${OUT_ROOT}/01_fastp" -type f \
      \( -name '*_merged.fastq.gz' -o -name '*_merged.fq.gz' \) 2>/dev/null \
      | sort \
      | while read -r merged_path; do
          sample_id="${merged_path##*/}"
          sample_id="${sample_id%_merged.fastq.gz}"
          sample_id="${sample_id%_merged.fq.gz}"
          [[ -n "${allow[$sample_id]:-}" ]] || continue
          printf "%s\t%s\n" "${sample_id}" "${merged_path}" >> "${PRINSEQ_INPUTS}"
        done
  fi

  # samples list (one SAMPLE per line)
  cut -f1 "${PRINSEQ_INPUTS}" | sort -u > "${PRINSEQ_LIST}" || :
  echo "[INFO] PRINSEQ inputs: $(wc -l < "${PRINSEQ_LIST}" | tr -d ' ') -> ${PRINSEQ_LIST}"
}

# ------------------------------- KRAKEN2 -------------------------------
# OUT_ROOT/00_logs/samples.for_kraken.txt
KRAKEN_SAMPLE_LIST="${OUT_ROOT}/00_logs/samples.for_kraken.txt"

build_kraken_sample_list() {
  : > "${KRAKEN_SAMPLE_LIST}"

  # Override: pass through verbatim (can be SAMPLE or absolute merged fastq path)
  if [[ -n "${OVERRIDE_LIST_KRAKEN:-}" && -f "${OVERRIDE_LIST_KRAKEN}" ]]; then
    awk 'NF>0 && $0 !~ /^#/' "${OVERRIDE_LIST_KRAKEN}" \
      | sed 's/[[:space:]]*$//' \
      | sort -u > "${KRAKEN_SAMPLE_LIST}"
    echo "[INFO] Kraken samples (override): $(wc -l < "${KRAKEN_SAMPLE_LIST}" | tr -d " ")"
    return
  fi

  # No override: derive SAMPLE IDs from PRIMARY (don’t scan disk)
  while read -r r; do
    [[ -z "$r" ]] && continue
    s="${r##*/}"
    s="${s%_R1_001.fastq.gz}"; s="${s%_R2_001.fastq.gz}"
    s="${s%_R1.fastq.gz}";     s="${s%_R2.fastq.gz}"
    s="${s%_1.fastq.gz}";      s="${s%_2.fastq.gz}"
    s="${s%_R1_001.fq.gz}";    s="${s%_R2_001.fq.gz}"
    s="${s%_R1.fq.gz}";        s="${s%_R2.fq.gz}"
    s="${s%_1.fq.gz}";         s="${s%_2.fq.gz}"
    [[ -n "$s" ]] && echo "$s"
  done < "${PRIMARY_LIST_PATH}" | sort -u > "${KRAKEN_SAMPLE_LIST}"

  echo "[INFO] Kraken samples: $(wc -l < "${KRAKEN_SAMPLE_LIST}" | tr -d " ") -> ${KRAKEN_SAMPLE_LIST}"
}

# ------------------------------- MAPPING (Bowtie2) -------------------------------
# We will build a list of SAMPLE IDs for mapping. The mapping script itself
# reads SAMPLE_LIST and uses SLURM_ARRAY_TASK_ID to pick the sample.

MAP_SAMPLE_LIST="${OUT_ROOT}/00_logs/samples.for_mapping.txt"

build_mapping_sample_list() {
  : > "${MAP_SAMPLE_LIST}"

  # Override: accept SAMPLE IDs or absolute *_GTDB_unclas.*q.gz; normalize to SAMPLE
  if [[ -n "${OVERRIDE_LIST_MAPPING:-}" && -f "${OVERRIDE_LIST_MAPPING}" ]]; then
    awk 'NF>0 && $0 !~ /^#/' "${OVERRIDE_LIST_MAPPING}" \
    | sed 's/[[:space:]]*$//' \
    | while IFS= read -r line; do
        if [[ "$line" == /* && "$line" == *_GTDB_unclas.*q.gz ]]; then
          b="${line##*/}"
          s="${b%_GTDB_unclas.fastq.gz}"
          s="${s%_GTDB_unclas.fq.gz}"
          printf '%s\n' "$s"
        else
          printf '%s\n' "$line"
        fi
      done \
    | sort -u > "${MAP_SAMPLE_LIST}"
    echo "[INFO] Mapping samples (override): $(wc -l < "${MAP_SAMPLE_LIST}" | tr -d " ")"
    return
  fi

  # No override: derive SAMPLE IDs from PRIMARY (don’t scan disk)
  while read -r r; do
    [[ -z "$r" ]] && continue
    s="${r##*/}"
    s="${s%_R1_001.fastq.gz}"; s="${s%_R2_001.fastq.gz}"
    s="${s%_R1.fastq.gz}";     s="${s%_R2.fastq.gz}"
    s="${s%_1.fastq.gz}";      s="${s%_2.fastq.gz}"
    s="${s%_R1_001.fq.gz}";    s="${s%_R2_001.fq.gz}"
    s="${s%_R1.fq.gz}";        s="${s%_R2.fq.gz}"
    s="${s%_1.fq.gz}";         s="${s%_2.fq.gz}"
    [[ -n "$s" ]] && echo "$s"
  done < "${PRIMARY_LIST_PATH}" | sort -u > "${MAP_SAMPLE_LIST}"

  echo "[INFO] Mapping samples: $(wc -l < "${MAP_SAMPLE_LIST}" | tr -d " ") -> ${MAP_SAMPLE_LIST}"
}

# ------------------------------- FILTERING / ngsLCA (+ bamdam) -------------------------------

FILTER_LIST="${OUT_ROOT}/00_logs/samples.for_filter.txt"

build_filter_sample_list() {
  : > "${FILTER_LIST}"

  # Override: keep lines verbatim (SAMPLE or absolute *.all.sorted.bam)
  if [[ -n "${OVERRIDE_LIST_FILTER:-}" && -f "${OVERRIDE_LIST_FILTER}" ]]; then
    awk 'NF>0 && $0 !~ /^#/' "${OVERRIDE_LIST_FILTER}" \
      | sed 's/[[:space:]]*$//' \
      | sort -u > "${FILTER_LIST}"
    echo "[INFO] Filtering items (override): $(wc -l < "${FILTER_LIST}" | tr -d " ")"
    return
  fi

  # No override: derive SAMPLE IDs from PRIMARY (don’t scan disk)
  while read -r r; do
    [[ -z "$r" ]] && continue
    s="${r##*/}"
    s="${s%_R1_001.fastq.gz}"; s="${s%_R2_001.fastq.gz}"
    s="${s%_R1.fastq.gz}";     s="${s%_R2.fastq.gz}"
    s="${s%_1.fastq.gz}";      s="${s%_2.fastq.gz}"
    s="${s%_R1_001.fq.gz}";    s="${s%_R2_001.fq.gz}"
    s="${s%_R1.fq.gz}";        s="${s%_R2.fq.gz}"
    s="${s%_1.fq.gz}";         s="${s%_2.fq.gz}"
    [[ -n "$s" ]] && echo "$s"
  done < "${PRIMARY_LIST_PATH}" | sort -u > "${FILTER_LIST}"

  echo "[INFO] Filtering samples: $(wc -l < "${FILTER_LIST}" | tr -d " ") -> ${FILTER_LIST}"
}

write_step_list() {
  local outfile="$1" force="$2" step="$3"
  local tmp="${outfile}.tmp"; : > "${tmp}"; cat > "${tmp}"
  if [[ "${force}" -eq 0 ]]; then
    case "${step}" in
      fastp)  awk -v out="${OUT_ROOT}/01_fastp"           '{s=$0; f=out"/"s"/"s"_merged.fastq.gz";                     cmd="test -f \""f"\""; if (system(cmd)!=0) print s}' "${tmp}" > "${outfile}" ;;
      sga)    awk -v out="${OUT_ROOT}/02_sga"             '{s=$0; f=out"/"s"/"s"_merged.dust.rmdup.fastq.gz";         cmd="test -f \""f"\""; if (system(cmd)!=0) print s}' "${tmp}" > "${outfile}" ;;
      kraken) awk -v out="${OUT_ROOT}/03_kraken_gtdb"     '{s=$0; f=out"/"s"_GTDB_unclas.fastq.gz";                   cmd="test -f \""f"\""; if (system(cmd)!=0) print s}' "${tmp}" > "${outfile}" ;;
      mapping)awk -v out="${OUT_ROOT}/04_mapping/''" '{s=$0; f=out"/"s".b2.k1000.all.sorted.bam";          cmd="test -f \""f"\""; if (system(cmd)!=0) print s}' "${tmp}" > "${outfile}" ;;
      filter) cp -f "${tmp}" "${outfile}" ;;
      *)      cp -f "${tmp}" "${outfile}" ;;
    esac
  else
    cp -f "${tmp}" "${outfile}"
  fi
  rm -f "${tmp}"
}

# -----------------------------------------------------------------------------
build_primary_list

# --- FASTP: discover pairs, build lists, submit with params ---
declare -A JID_FASTP
FASTP_PAIRS="${OUT_ROOT}/00_logs/fastp.pairs.tsv"      # SAMPLE \t R1 \t R2
FASTP_TODO="${OUT_ROOT}/00_logs/samples_for_fastp.txt" # SAMPLE (for dependency visibility)

pair_fastq_from_primary() {
  : > "${FASTP_PAIRS}"
  while read -r f; do
    [[ -z "$f" ]] && continue
    # only seed off R1/_1; skip R2 lines
    case "$f" in
      *R1*.fastq.gz|*R1*.fq.gz)
        r1="$f"; r2="${f/R1/R2}"
        ;;
      *_1*.fastq.gz|*_1*.fq.gz)
        r1="$f"; r2="${f/_1/_2}"
        ;;
      *) continue ;;
    esac
    [[ -f "$r2" ]] || { echo "[WARN] No mate for $r1"; continue; }
    # sample name = basename without the R1 token + tail
    s="$(basename "$r1")"
    s="${s%_R1_001.fastq.gz}"; s="${s%_R1.fastq.gz}"
    s="${s%_R1_001.fq.gz}";   s="${s%_R1.fq.gz}"
    s="${s%_1.fastq.gz}";     s="${s%_1.fq.gz}"
    printf "%s\t%s\t%s\n" "$s" "$r1" "$r2" >> "${FASTP_PAIRS}"
  done < "${PRIMARY_LIST_PATH}"
  # build todo list (skip if merged exists unless forcing)
  : > "${FASTP_TODO}"
  while IFS=$'\t' read -r s r1 r2; do
    [[ -z "$s" ]] && continue
    out="${OUT_ROOT}/01_fastp/${s}/${s}_merged.fastq.gz"
    if [[ ${FORCE_FASTP:-0} -eq 1 || ! -f "$out" ]]; then
      echo "$s" >> "${FASTP_TODO}"
    fi
  done < "${FASTP_PAIRS}"
}
if [[ ${ENABLE_PREPROCESS:-1} -eq 1 && ${ENABLE_FASTP:-1} -eq 1 ]]; then
  echo "[INFO] FASTP: pairing from ${PRIMARY_LIST_PATH}"
  pair_fastq_from_primary
  n_fastp=$(wc -l < "${FASTP_TODO}" | tr -d ' ' || echo 0)
  echo "[INFO] FASTP todo: ${n_fastp}"
  # index pairs once for quick lookup
  declare -A R1_BY_S R2_BY_S
  while IFS=$'\t' read -r s r1 r2; do
    [[ -n "$s" ]] && { R1_BY_S["$s"]="$r1"; R2_BY_S["$s"]="$r2"; }
  done < "${FASTP_PAIRS}"
  while read -r sname; do
    [[ -z "$sname" ]] && continue
    r1="${R1_BY_S[$sname]}"; r2="${R2_BY_S[$sname]}"
    [[ -f "$r1" && -f "$r2" ]] || { echo "[WARN] Missing R1/R2 for $sname"; continue; }
    outdir="${OUT_ROOT}/01_fastp/${sname}"; mkdir -p "$outdir"
    sbatch_opts FASTP
    jid=$(sbatch "${SBATCH_BUILT_OPTS[@]}" \
      --job-name="${FASTP_SBATCH_JOB_NAME}_${sname}" \
      --output="${LOG_ROOT}/out/fastp.%x.%j.out" \
      --error="${LOG_ROOT}/error/fastp.%x.%j.err" \
      --export=ALL,\
SAMPLE="${sname}",\
R1_IN="${r1}",\
R2_IN="${r2}",\
OUTPUT="${outdir}",\
FASTP_MODULE="${FASTP_MODULE}",\
FASTP_OVERLAP_LEN_REQUIRE="${FASTP_OVERLAP_LEN_REQUIRE}",\
FASTP_MIN_LENGTH="${FASTP_MIN_LENGTH}" \
      "${SCRIPTS_DIR}/01_fastp.sh" | awk '{print $4}')
    JID_FASTP["$sname"]="$jid"
    echo "[FASTP] ${sname} -> ${jid}"
  done < "${FASTP_TODO}"
else
  echo "[INFO] FASTP disabled."
fi

# -----------------------------------------------------------------------------
# SGA
# -----------------------------------------------------------------------------

build_sga_inputs

# --- SGA submission ---
declare -A JID_SGA
if [[ ${ENABLE_PREPROCESS:-1} -eq 1 && ${ENABLE_SGA:-1} -eq 1 ]]; then
  echo "[INFO] Submitting SGA…"
  sbatch_opts SGA
  while IFS=$'\t' read -r s merged; do
    [[ -z "$s" || -z "$merged" ]] && continue
    outdir="${OUT_ROOT}/02_sga/${s}"; mkdir -p "${outdir}"

    dep=()
    if [[ -n "${JID_FASTP[$s]:-}" ]]; then
      dep=( --dependency="afterok:${JID_FASTP[$s]}" )
    fi

    jid=$(sbatch "${dep[@]}" "${SBATCH_BUILT_OPTS[@]}" \
      --job-name="${SGA_SBATCH_JOB_NAME}_${s}" \
      --output="${LOG_ROOT}/out/fastp.%x.%j.out" \
      --error="${LOG_ROOT}/error/fastp.%x.%j.err" \
      --export=ALL,\
SAMPLE="${s}",\
INPUT_MERGED="${merged}",\
OUTPUT_DIR="${outdir}",\
CONDA_INIT="${CONDA_INIT}",\
CONDA_ENV_SGA="${CONDA_ENV_SGA}",\
SGA_DUST_THRESHOLD="${SGA_DUST_THRESHOLD}" \
      "${SCRIPTS_DIR}/02_sga.sh" | awk '{print $4}')
    JID_SGA["$s"]="$jid"
    echo "[SGA] ${s} -> ${jid}"
  done < "${SGA_INPUTS}"
else
  echo "[INFO] SGA disabled."
fi


# ---------------------------------------------------------------------------------------------
# PRINSEQ (alternative to SGA)
# ---------------------------------------------------------------------------------------------

typeset -p JOB_PRINSEQ &>/dev/null || JOB_PRINSEQ=""
if [[ ${ENABLE_PREPROCESS:-1} -eq 1 && ${ENABLE_PRINSEQ:-0} -eq 1 ]]; then
  build_prinseq_inputs
  n_prinseq=$(wc -l < "${PRINSEQ_LIST}" | tr -d ' ' || echo 0)
  if [[ "${n_prinseq}" -eq 0 ]]; then
    echo "[INFO] PRINSEQ: nothing to do."
  else
    # Depend on FASTP per this invocation (if we launched it)
    # depend on any FASTP jobs launched in this invocation
    deps=()
    if ((${#JID_FASTP[@]} > 0)); then
      dep_ids=$(printf "%s," "${JID_FASTP[@]}" | sed 's/,$//')
      [[ -n "${dep_ids}" ]] && deps=( --dependency="afterok:${dep_ids}" )
    fi

    sbatch_opts PRINSEQ
    require_file "${SCRIPTS_DIR}/03_prinseq.sh"
    mkdir -p "${OUT_ROOT}/02_prinseq"

    JOB_PRINSEQ=$(sbatch "${deps[@]}" "${SBATCH_BUILT_OPTS[@]}" \
      --job-name="${PRINSEQ_SBATCH_JOB_NAME:-prinseq}_arr" \
      --array="0-$((n_prinseq-1))" \
      --output="${LOG_ROOT}/out/fastp.%x.%j.out" \
      --error="${LOG_ROOT}/error/fastp.%x.%j.err" \
      --export=ALL,\
SAMPLE_LIST="${PRINSEQ_LIST}",\
PRINSEQ_INPUTS_TSV="${PRINSEQ_INPUTS}",\
FASTP_DIR="${OUT_ROOT}/01_fastp",\
OUT_DIR="${OUT_ROOT}/02_prinseq",\
PRINSEQ_LITE="${PRINSEQ_LITE}",\
PRINSEQ_COMPLEXITY_METHOD="${PRINSEQ_COMPLEXITY_METHOD}",\
PRINSEQ_COMPLEXITY_THRESHOLD="${PRINSEQ_COMPLEXITY_THRESHOLD}",\
PRINSEQ_MIN_LEN="${PRINSEQ_MIN_LEN}",\
PRINSEQ_DEREP="${PRINSEQ_DEREP}" \
      "${SCRIPTS_DIR}/03_prinseq.sh" | awk '{print $4}')
    echo "[PRINSEQ] array -> ${JOB_PRINSEQ} (0-$((n_prinseq-1)))"
  fi
fi

# -----------------------------------------------------------------------------
# Kraken2 GTDB (single job)
# At that step we collect array job ID of all SGA or prinseq and run one Kraken2 GTDB job
# So the logic is to wait that all SGA/prinseq are finished before to run Kraken2
# -----------------------------------------------------------------------------
# Ensure JID_SGA exists even if SGA step didn’t run in this invocation
if ! declare -p JID_SGA &>/dev/null; then
  declare -A JID_SGA=()
fi

if [[ ${ENABLE_KRAKEN_GTDB:-1} -eq 1 ]]; then
  # Build sample list (override or discovered)
  if [[ -n "${OVERRIDE_LIST_KRAKEN:-}" && -f "${OVERRIDE_LIST_KRAKEN}" ]]; then
    awk 'NF>0 {print $1}' "${OVERRIDE_LIST_KRAKEN}" | sort -u > "${KRAKEN_SAMPLE_LIST}"
    echo "[INFO] Kraken samples (override): $(wc -l < "${KRAKEN_SAMPLE_LIST}" | tr -d ' ')"
  else
    build_kraken_sample_list
  fi

  # Choose input source + dependency
  deps=()
  if [[ ${ENABLE_PRINSEQ:-0} -eq 1 ]]; then
    KRAKEN_INPUT_DIR="${OUT_ROOT}/02_prinseq"
    KRAKEN_INPUT_MODE="prinseq"
    [[ -n "${JOB_PRINSEQ:-}" ]] && deps=( --dependency="afterok:${JOB_PRINSEQ}" )
  else
    KRAKEN_INPUT_DIR="${OUT_ROOT}/02_sga"
    KRAKEN_INPUT_MODE="sga"
    # Accept either a scalar JOB_SGA or your existing associative array JID_SGA[*]
    if [[ -n "${JOB_SGA:-}" ]]; then
      deps=( --dependency="afterok:${JOB_SGA}" )
    elif ((${#JID_SGA[@]} > 0)); then
      dep_ids=$(printf "%s:" "${JID_SGA[@]}" | sed 's/:$//')
      [[ -n "${dep_ids}" ]] && deps=( --dependency="afterok:${dep_ids}" )
    fi
  fi

  sbatch_opts KRAKEN
  mkdir -p "${OUT_ROOT}/03_kraken_gtdb"
  require_file "${SCRIPTS_DIR}/04_Kraken_gtdb.sh"

  JID_KRAKEN=$(
    sbatch "${deps[@]}" "${SBATCH_BUILT_OPTS[@]}" \
      --job-name="${KRAKEN_SBATCH_JOB_NAME}_batch" \
      --output="${LOG_ROOT}/out/kraken.%x.%j.out" \
      --error="${LOG_ROOT}/error/kraken.%x.%j.err" \
      --export=ALL,\
SAMPLE_LIST="${KRAKEN_SAMPLE_LIST}",\
INPUT_DIR="${KRAKEN_INPUT_DIR}",\
INPUT_MODE="${KRAKEN_INPUT_MODE}",\
OUTPUT_DIR="${OUT_ROOT}/03_kraken_gtdb",\
KRAKEN2_MODULE="${KRAKEN2_MODULE}",\
PDC_MODULE="${PDC_MODULE}",\
GTDB_SRC="${GTDB_SRC}",\
THREADS="${KRAKEN_SBATCH_CPUS}" \
      "${SCRIPTS_DIR}/04_Kraken_gtdb.sh" | awk '{print $4}'
  )
  echo "[KRAKEN] batch -> ${JID_KRAKEN}"
else
  echo "[INFO] Kraken disabled."
  JID_KRAKEN=""
fi

# -----------------------------------------------------------------------------
# Mapping (array)
# -----------------------------------------------------------------------------

if [[ ${ENABLE_MAPPING:-1} -eq 1 ]]; then
  build_mapping_sample_list

  n_map=$(wc -l < "${MAP_SAMPLE_LIST}" | tr -d ' ' || echo 0)
  if [[ "${n_map}" -eq 0 ]]; then
    echo "[INFO] Mapping: nothing to do."
  else
    # Depend on Kraken batch if we launched it here
    deps=()
    if [[ -n "${JID_KRAKEN:-}" ]]; then
      deps=( --dependency="afterok:${JID_KRAKEN}" )
    fi

    sbatch_opts MAP
    require_file "${SCRIPTS_DIR}/05_mapping_bowtie2.sh"

    jid_map=$(sbatch "${deps[@]}" "${SBATCH_BUILT_OPTS[@]}" \
      --job-name="${MAP_SBATCH_JOB_NAME:-bowtie2}_arr" \
      --array="0-$((n_map-1))" \
      --output="${LOG_ROOT}/out/fastp.%x.%j.out" \
      --error="${LOG_ROOT}/error/fastp.%x.%j.err" \
      --export=ALL,\
SAMPLE_LIST="${MAP_SAMPLE_LIST}",\
INPUT_DIR="${OUT_ROOT}/03_kraken_gtdb",\
OUTPUT_DIR="${OUT_ROOT}/04_mapping",\
TMPDIR="${TMP_ROOT}",\
PHYLONORWAY="${PHYLONORWAY}",\
HEADER="${PHYLONORWAY_HEADER}",\
PLASTID="${PLASTID}",\
MITO="${MITO}",\
MAM_BIRD_FISH="${MAM_BIRD_FISH}",\
NAMES="${NAMES}",\
NODES="${NODES}",\
ACC2TAX="${ACC2TAX}",\
BOWTIE2_MODULE="${BOWTIE2_MODULE}",\
SAMTOOLS_MODULE="${SAMTOOLS_MODULE}" \
      "${SCRIPTS_DIR}/05_mapping_bowtie2.sh" | awk '{print $4}')
    JOB_MAP="${jid_map}"     # <-- set ONLY when we submit mapping now
    echo "[MAPPING] array -> ${jid_map} (0-$((n_map-1)))"
  fi
else
  echo "[INFO] Mapping disabled."
fi

# ----------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Filtering + ngsLCA (single job)
# -----------------------------------------------------------------------------
: "${JOB_MAP:=}"   # mapping array job id if mapping ran in this invocation

FILTER_LIST="${OUT_ROOT}/00_logs/samples.for_filter.txt"
if [[ ${ENABLE_FILTERING:-1} -eq 1 ]]; then
  build_filter_sample_list
  n_filter=$(wc -l < "${FILTER_LIST}" | tr -d ' ' || echo 0)
  if [[ "${n_filter}" -eq 0 ]]; then
    echo "[INFO] Filtering: nothing to do."
  else
    # Depend on the mapping array if present
    deps=()
    [[ -n "${JOB_MAP}" ]] && deps=( --dependency="afterok:${JOB_MAP}" )

    sbatch_opts FILTER
    require_file "${SCRIPTS_DIR}/06_filtering.sh"
    mkdir -p "${OUT_ROOT}/05_filtering" "${OUT_ROOT}/05_filtering/bamdam"

    JOB_FILTER=$(sbatch "${deps[@]}" "${SBATCH_BUILT_OPTS[@]}" \
      --job-name="${FILTER_SBATCH_JOB_NAME:-filtering_ngslca}_arr" \
      --array="0-$((n_filter-1))" \
      --output="${LOG_ROOT}/out/fastp.%x.%j.out" \
      --error="${LOG_ROOT}/error/fastp.%x.%j.err" \
      --export=ALL,\
SAMPLE_LIST="${FILTER_LIST}",\
MAP_DIR="${OUT_ROOT}/04_mapping",\
OUT_DIR="${OUT_ROOT}/05_filtering",\
ENABLE_FILTERBAM="${ENABLE_FILTERBAM}",\
ENABLE_NGSLCA="${ENABLE_NGSLCA}",\
ENABLE_BAMDAM="${ENABLE_BAMDAM}",\
CONDA_ENV_FILTERBAM="${CONDA_ENV_FILTERBAM}",\
CONDA_ENV_ngsLCA="${CONDA_ENV_ngsLCA}",\
NAMES="${NAMES}",\
NODES="${NODES}",\
ACC2TAX="${ACC2TAX}",\
BAMDAM_PYTHON_MODULE="${BAMDAM_PYTHON_MODULE}",\
BAMDAM_VENV="${BAMDAM_VENV}",\
KRONATOOLS_MODULE="${KRONATOOLS_MODULE}",\
THREADS="${FILTER_SBATCH_CPUS}",\
BAMDAM_STRANDED="${BAMDAM_STRANDED}",\
BAMDAM_MINREADS="${BAMDAM_MINREADS}",\
BAMDAM_MAXDAMAGE="${BAMDAM_MAXDAMAGE}" \
      "${SCRIPTS_DIR}/06_filtering.sh" | awk '{print $4}')
    echo "[FILTER] array -> ${JOB_FILTER} (0-$((n_filter-1)))"
  fi
else
  echo "[INFO] Filtering/ngsLCA disabled."
fi
# -------------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# MMSeqs2 (optional stub)
# -----------------------------------------------------------------------------
if [[ ${ENABLE_MMSEQS2:-0} -eq 1 ]]; then
  echo "[INFO] MMSeqs2 is enabled; add commands in "${SCRIPTS_DIR}/07_mmseqs2.sh"."
fi

# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------
last_dep=""
for jid in "${job_filter:-}" "${job_map:-}" "${job_kraken:-}"; do
  if [[ -n "${jid}" ]]; then last_dep="${jid}"; break; fi
done

mkdir -p "${OUT_ROOT}/99_metrics"
sbatch_opts METRICS
if [[ -n "${last_dep}" ]]; then
  sbatch --dependency=afterok:${last_dep} "${SBATCH_BUILT_OPTS[@]}" \
    --export=ALL,ROOT="${OUT_ROOT}",MAP_RUN="${MAP_RUN}" \
    "${SCRIPTS_DIR}/99_metrics_to_tsv.sh" >/dev/null
else
  sbatch "${SBATCH_BUILT_OPTS[@]}" \
    --export=ALL,ROOT="${OUT_ROOT}",MAP_RUN="${MAP_RUN}" \
    "${SCRIPTS_DIR}/99_metrics_to_tsv.sh" >/dev/null
fi

echo "[DONE] Submission complete. Check ${LOG_ROOT}/out and ${LOG_ROOT}/error."
