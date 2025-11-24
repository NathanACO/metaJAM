#!/bin/bash -l
#SBATCH -A naiss2025-5-172
#SBATCH -p shared
#SBATCH --cpus-per-task=4
#SBATCH --time 02:00:00
#SBATCH --job-name=metrics
#SBATCH -e /cfs/klemming/projects/supr/naiss2025-23-301/output/error/%x-%j.error
#SBATCH -o /cfs/klemming/projects/supr/naiss2025-23-301/output/out/%x-%j.out

set -euo pipefail

ROOT="${ROOT}"
MAP_RUN="${MAP_RUN}"
OUT="${ROOT}/99_metrics"
mkdir -p "${OUT}"

echo -e "sample\tfastp_merged_reads\tpostSGA_reads\tGTDB_unclas_reads\tbam_total_reads" > "${OUT}/metrics.tsv"

# Count reads without changing your upstream commands:
for sdir in "${ROOT}/01_fastp"/*; do
  [[ -d "$sdir" ]] || continue
  sample=$(basename "$sdir")
  fastp_mgz="${ROOT}/01_fastp/${sample}/${sample}_merged.fastq.gz"
  sga_mgz="${ROOT}/02_sga/${sample}/${sample}_merged.dust.rmdup.fastq.gz"
  gtdb_unclas="${ROOT}/03_kraken_gtdb/${sample}_GTDB_unclas.fastq.gz"
  bam="${ROOT}/04_mapping/${MAP_RUN}/${sample}.b2.k1000.all.sorted.bam"

  # counts (fastq: lines/4; bam: use samtools idxstats if index exists, else flagstat reads)
  f_count=$( (zcat "$fastp_mgz" 2>/dev/null || gzip -dc "$fastp_mgz") | awk 'END{print NR/4}' )
  s_count=$( (zcat "$sga_mgz" 2>/dev/null || gzip -dc "$sga_mgz") | awk 'END{print NR/4}' )
  u_count=$( (zcat "$gtdb_unclas" 2>/dev/null || gzip -dc "$gtdb_unclas") | awk 'END{print NR/4}' )
  b_count=$(samtools view -c "$bam" 2>/dev/null || echo 0)

  echo -e "${sample}\t${f_count}\t${s_count}\t${u_count}\t${b_count}" >> "${OUT}/metrics.tsv"
done

echo "Wrote ${OUT}/metrics.tsv  (Excel can open TSV directly)."