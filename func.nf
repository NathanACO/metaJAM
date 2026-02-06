#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process FASTP {

    label 'small_memory'

    conda 'bioconda::fastp'

    input:
    tuple val(ID), path(r1), path(r2), val(FASTP_OVERLAP_LEN_REQUIRE), val(FASTP_MIN_LENGTH)

    output:
    tuple val(ID), path("*_merged.fastq.gz")

    publishDir "${params.OUTPUT_Dir}/01_fastp", mode: "copy" 

    script:
    """

        fastp -i "${r1}" -I "${r2}" -p -c --merge \
        --overlap_len_require "${FASTP_OVERLAP_LEN_REQUIRE}" \
        --overlap_diff_limit 1 \
        --merged_out="${ID}_merged.fastq.gz" \
        --include_unmerged \
        --trim_poly_x \
        -h "${ID}.html" -R "${ID} QC Report" \
        -w 10 \
        -l "${FASTP_MIN_LENGTH}"

    """
}


process SGA {
    label 'small_memory'
	conda ' bioconda::sga'
	input:
    	tuple val(ID), path(merged_reads), val(SGA_DUST_THRESHOLD)

    output:
		tuple val(ID), path("*_merged.dust.rmdup.fastq.gz"), emit: "rm_low_complexity"
		path("*_merged.dust.fastq.gz")

    publishDir "${params.OUTPUT_Dir}/02_sga", mode: "copy" 

    script:
    """
		# Preprocess + index + filter
		sga preprocess "${merged_reads}" \
		-o "${ID}_merged.dust.fastq.gz" \
		-m 30 \
		--dust \
		--dust-threshold="${SGA_DUST_THRESHOLD}" \
		--no-primer-check \
		--min-length=30 | tee -a "sga_${ID}.log"

		sga index "${ID}_merged.dust.fastq.gz" \
		--algorithm=ropebwt \
		--threads="${task.cpus}" \
		-p "${ID}_merged.dust"

		sga filter --homopolymer-check \
		--no-kmer-check \
		--threads "${task.cpus}" \
		-p "${ID}_merged.dust" \
		"${ID}_merged.dust.fastq.gz" \
		-o "${ID}_merged.dust.rmdup.fastq.gz"

    """
}


process PRINSEQ {
    label 'small_memory'
	conda 'bioconda::prinseq'
	input:
    tuple val(ID), path(merged_reads), val(PRINSEQ_COMPLEXITY_METHOD), val(PRINSEQ_COMPLEXITY_THRESHOLD), val(PRINSEQ_MIN_LEN), val(PRINSEQ_DEREP)

    output:
    tuple val(ID), path("*_merged.complexity_filtered.duplicatesremoved*"), emit: "rm_low_complexity"    
    tuple path("*_merged.complexity_filtered*"), path("*_merged.low_complexity*")

    publishDir "${params.OUTPUT_Dir}/02_prinseq", mode: "copy" 

    script:
    """
		# File stems (NO per-sample folder)
		GOOD="${ID}_merged.complexity_filtered"
		BAD="${ID}_merged.low_complexity"
		DEREP_OUT="${ID}_merged.complexity_filtered.duplicatesremoved"

		# 1) Low complexity filter
		zcat "${merged_reads}" | prinseq-lite.pl \
		-fastq stdin \
		-out_good "\${GOOD}" \
		-out_bad "\${BAD}" \
		-lc_method "${PRINSEQ_COMPLEXITY_METHOD}" \
		-lc_threshold "${PRINSEQ_COMPLEXITY_THRESHOLD}" \
		-min_len "${PRINSEQ_MIN_LEN}" \
		-line_width 0

		gzip -f "\${BAD}.fastq"

		# 3) Remove exact fwd/rev duplicates
		prinseq-lite.pl \
		-fastq "\${GOOD}.fastq" \
		-out_good "\${DEREP_OUT}" \
		-out_bad null \
		-min_len "${PRINSEQ_MIN_LEN}" \
		-derep "${PRINSEQ_DEREP}" \
		-line_width 0

		gzip -f "\${DEREP_OUT}.fastq"
		"""
}

process KRAKEN2 {
    label 'small_memory' // for test
    conda 'bioconda::kraken2'

    publishDir "${params.OUTPUT_Dir}/03_kraken2_filter", mode: "copy" 

    input:    
        tuple val(ID), path(reads), path(DB)

    output:
        tuple val(ID), path("*_unclas.fastq.gz"), emit: not_microbe
        path("*_report.txt")
        path("*_output.txt*")
        
    script:
    """
        DB_LABEL=\$(basename $DB)
        kraken2 --db "${DB}" \
            --report "${ID}_\${DB_LABEL}_report.txt" \
            --report-minimizer-data \
            --gzip-compressed \
            --threads ${task.cpus} \
            --output "${ID}_\${DB_LABEL}_output.txt" \
            --classified-out "${ID}_\${DB_LABEL}_clas.fastq" \
            --unclassified-out "${ID}_\${DB_LABEL}_unclas.fastq" \
            --memory-mapping "${reads}"

        gzip *.fastq
        gzip "${ID}_\${DB_LABEL}_output.txt"
    """
}


process BOWTIE2 {
    label 'small_memory' // for test
    conda './envs/bowtie2.yml' // to check
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)
              //path(idx1), path(idx2),path(idx3), path(idx4),path(idx5), path(idx6) // unlimited number of DBs

output:
    tuple val(ID), path("*.bam")

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
        bowtie2 --very-sensitive -p "${task.cpus}" -k $n_allow_multimapper -x \$(basename $idx) -U $reads | \
            samtools view -@ "${task.cpus}" -Sb -q 0 -F 4 - > $ID.\$(basename $idx).bam
    """
}

process GET_ACC2TAXID {
    label 'little_memory'

    conda ''

    input:
    tuple val(idx_ID), path(idxs)

    output:
    path "*.acc2taxid"

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"

    script:
    """
        if acc2taxid_exists == "false"; do
            bowtie2-inspect $idx_ID | grep ">" | cut -f 1,2,3 -d" "| sed 's/>//' > $idx_ID.contigs
        fi

        generate_acc2taxid.sh $idx_ID.contigs ${task.cpus}
    """


}

// merge bam of same sample mapped to different databases
process MERGE_BAM {
    label 'little_memory'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(bams)

    output:
        tuple val(ID), path("*_merged.sorted.bam")
    
    publishDir "${params.OUTPUT_Dir}/05_merged_bam", mode: "copy"
        
    script:
    """
        #filtering out unmapped reads in case it's not done for input bam
        for bam1 in *.bam; do
            [[ "\$bam1" == *mapped.bam ]] && continue
            samtools view -@ "${task.cpus}" -b -F 0x4 "\$bam1" -o "\$(basename \$bam1 .bam).mapped.bam"
        done

        samtools merge ${ID}.merged.bam *.mapped.bam

        samtools quickcheck ${ID}.merged.bam || {
            echo "ERROR: Merging is not successful: ${ID}.merged.bam" >&2
            exit 1
        }

        samtools sort -@ "${task.cpus}" -n   -o ${ID}_merged.sorted.bam   ${ID}.merged.bam
    """
}



// merge bam of same sample mapped to different databases
process MASK_REGIONS {
    label 'little_memory'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(bam), path(bed)

    output:
        tuple val(ID), path("*.mic_masked.bam")
    
    publishDir "${params.OUTPUT_Dir}/06_masked_bam", mode: "copy"
        
    script:
    """
        samtools view -b -h -L "$bed" -U "\$(basename $bam .bam).mic_masked.bam" -o /dev/null "$bam" 
    """
}

process FILTERBAM {

    label 'little_memory'
    conda './envs/filterBAM.yml'
    input:    
        tuple val(ID), path(bam)
            
    output:
        tuple val(ID), path("*.b2.k1000.all.filtered.bam")

    publishDir "${params.OUTPUT_Dir}/filterbam", mode: "copy"
        
    script:
    """
        filterBAM filter \
            --bam "${bam}" \
            --bam-filtered "${ID}.b2.k1000.all.filtered.bam" \
            --stats "${ID}.b2.k1000.all.stats.tsv.gz" \
            --stats-filtered "${ID}.b2.k1000.all.stats-filtered.tsv.gz" \
            --threads "${task.cpus}" \
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
    """
}

process NGSLCA {

    label 'little_memory'
    conda './envs/ngsLCA2.yml'
    input:    
        tuple val(ID), path(bam), path(NAMES), path(NODES), path(ACC2TAX)
            
    output:
        tuple val(ID), path("*.lca")

    publishDir "${params.OUTPUT_Dir}/07_ngslca", mode: "copy"
        
    script:
    """
        ngsLCA \
            -simscorelow 0.95 \
            -simscorehigh 1 \
            -names "${NAMES}" \
            -nodes "${NODES}" \
            -acc2tax "${ACC2TAX}" \
            -bam "${bam}" \
            -outnames "\$(basename $bam .bam)" \
            -fix_ncbi 0
    """
}

process BAMDAM {
    conda './envs/bamdam.yml'

    publishDir "${params.OUTPUT_Dir}/08_bamdam", mode: "copy"

    cpus { 25 * task.attempt }
    memory { 20.GB * task.attempt }
    time { 2.hour }

    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' } // in case there is no reads
    maxRetries 3

    input:    
        tuple val(ID), path(bam), path(lca), 
        val(STRANDED), 
        val(MINREADS), 
        val(MAXDAMAGE), 
        val(TOP_GENUS)
            
    output:
        tuple val(ID), path("*.small.bam"), path("*.small.lca"), path("*.tsv"), emit: lca
        path("*.subs.txt")
        tuple val(ID), path("*.xml"), emit: xml
        
    publishDir "${params.OUTPUT_Dir}/bamdam", mode: "copy"

    script:
    """
        # shrink
        bamdam shrink \
        --in_bam "${bam}" \
        --in_lca "${lca}" \
        --out_bam "${ID}.small.bam" \
        --out_lca "${ID}.small.lca" \
        --stranded "${STRANDED}" \
        --show_progress

    if [ ! -s "${ID}.small.lca" ]; then
        echo "bamdam shrink failed"
        exit 1
    fi

    if [ "\$(wc -l < "${ID}.small.lca")" -eq 0 ]; then
        echo "no reads remained"
        exit 1
    fi

        # compute
        bamdam compute \
        --in_bam "${ID}.small.bam" \
        --in_lca "${ID}.small.lca" \
        --out_tsv "${ID}.tsv" \
        --out_subs "${ID}.subs.txt" \
        --stranded "${STRANDED}" \
        --show_progress

        # krona + HTML (per-sample)
        bamdam krona \
        --in_tsv "${ID}.tsv" \
        --out_xml "${ID}.xml" \
        --minreads "${MINREADS}" \
        --maxdamage "${MAXDAMAGE}"

        while read -r count taxid genus; do
        bamdam plotdamage \
        --in_subs "${ID}.subs.txt" \
        --tax "\$taxid" \
        --outplot "damageplot_\${genus}_\${taxid}_${ID}.png" \
        || echo "Warning: bamdam failed for taxid \$taxid (continuing)"
        done < <(
        awk '{for(i=2;i<=NF;i++){n=split(\$i,a,":");if(n>=3 && a[3]=="genus"){id=a[1];counts[id]++;name[id]=a[2];break}}}
        END{for(id in counts)printf "%d\t%s\t%s\n",counts[id],id,name[id]}' \
        "${ID}.small.lca" \
        | sort -nrk1,1 \
        | head -n "$TOP_GENUS"
        )
  """
}


process KRONATOOLS {
    label 'small_memory'
    conda 'bioconda::krona'
    publishDir "${params.OUTPUT_Dir}/09_kronaTools", mode: "copy"
    input:
        tuple val(ID), path(xml)
            
    output:
        path("*.html")

    script:
    """
        ktImportXML -o "${ID}.html" "$xml"
    """
}



process MMSEQ2 {
    conda 'envs/mmseq2.yml'

    label 'small_memory'

    // setting for NT database
    //cpus 256
    //memory 880.GB
    //time 10.hour
    //queue 'memory'

    
    input:    
        tuple val(ID), path(bam), path(lca), path(lca_tsv),
        val(MIN_DMG),
        val(TOP_GENERA),
        val(MAX_READS),
        val(MIN_READS),
        val(SPACED_KMER_MODE),
        val(S),
        val(MAX_EVALUE),
        val(MIN_LENGTH),
        val(MIN_SEQID),
        val(MAX_SEQS),
        val(MIN_QUERY_COV),
        val(SPLIT_MEM_LIMIT),
        val(AMBIG_FRAC),
        val(DATABASE_NAME),
        val(MMSEQS2_DB),
        val(SEED),
        val(MIN_BITS),
        path(TAXADB_SQLITE),
        path(GENERA_FILE)

    output:
        tuple val(ID),
        path("*.mmseqs_queries.fa"),
        path("*.mmseqs_expected.tsv"),
        path("*.besthit.assigned.tsv"),
        path("*.bamdam_mmseqs.evaluation.summary.tsv")


    publishDir "${params.OUTPUT_Dir}/bamdam", mode: "copy"

    script:
    def seedArg = SEED ? "--seed ${SEED}" : ""
    def genusArg = GENERA_FILE.name != 'NO_FILE' ? "--genera-file ${GENERA_FILE}" : ''

    """
    # -----------------------------------------------------------------------------
    # 1) Build combined FASTA + expected map from LCA (top10 genus / kingdom)
    # -----------------------------------------------------------------------------
    filter_lca_top10_subset_reads.py \
    --lca "${lca}" \
    --bamdam-tsv "${lca_tsv}" \
    --min-dmg "${MIN_DMG}" --top-genera "${TOP_GENERA}" \
    --max-reads "${MAX_READS}" \
    --min-reads "${MIN_READS}" \
    ${seedArg} \
    ${genusArg} \
    --out-fasta "${ID}.mmseqs_queries.fa" \
    --out-map   "${ID}.mmseqs_expected.tsv"

    # If nothing was written (all taxa skipped), stop gracefully
    if [[ ! -s "${ID}.mmseqs_queries.fa" ]]; then
    echo "[INFO] No queries generated for ${ID} (likely all taxa < min reads)."
    exit 0
    fi

    # -----------------------------------------------------------------------------
    # 2) MMseqs2 search
    # -----------------------------------------------------------------------------
    mmseqs createdb "${ID}.mmseqs_queries.fa" "${DATABASE_NAME}"

    # if [[ -f "${DATABASE_NAME}.resultDB.dbtype" ]]; then
    # echo "[INFO] MMseqs search output exists, reusing: ${DATABASE_NAME}.resultDB"
    # else
    mmseqs search \
        "${DATABASE_NAME}" \
        "${MMSEQS2_DB}" \
        "${DATABASE_NAME}.resultDB" \
        "./" \
        --threads "${task.cpus}" \
        -a \
        --search-type 3 \
        --spaced-kmer-mode "${SPACED_KMER_MODE}" \
        -s "${S}" \
        -e "${MAX_EVALUE}" \
        --min-length "${MIN_LENGTH}" \
        --min-seq-id "${MIN_SEQID}" \
        --max-seqs "${MAX_SEQS}" \
        --cov-mode 2 -c "${MIN_QUERY_COV}" \
        --split-memory-limit "${SPLIT_MEM_LIMIT}"
    # fi

    # if [[ -s "${DATABASE_NAME}.resultDB.m8" ]]; then
    # echo "[INFO] convertalis output exists, reusing: ${DATABASE_NAME}.resultDB.m8"
    # else
    mmseqs convertalis \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,tlen,qlen,tcov,qcov,qaln,taln" \
        "${DATABASE_NAME}" \
        "${MMSEQS2_DB}" \
        "${DATABASE_NAME}.resultDB" \
        "${DATABASE_NAME}.resultDB.m8" \
        --format-mode 4 \
        --search-type 3
    # fi

    ## Bits filtering ##
    awk -F'\t' "NR==1 || \$12 > ${MIN_BITS}" "${DATABASE_NAME}.resultDB.m8" > "${DATABASE_NAME}.resultDB_bits.m8"

    # -----------------------------------------------------------------------------
    # 3) Add taxonomic info to MMseqs2 m8 (optional; requires TAXADB)
    # -----------------------------------------------------------------------------
    add_taxid_info.py \
        -blst Mmseqs2 \
        -b "${DATABASE_NAME}.resultDB_bits.m8" \
        -d "${TAXADB_SQLITE}" \
        -o "${ID}.${DATABASE_NAME}.resultDB_bits_nucl_taxid.m8" \
        -del tab
    # deactivate
    # else
    # echo "[INFO] TAXADB not configured; skipping taxid augmentation."
    # fi

    # rm -f "nt_query.resultDB."[0-9]* # Remove intermediate files

    # -----------------------------------------------------------------------------
    # Post-processing: pick one best hit per query + compare taxonomy vs expected
    # -----------------------------------------------------------------------------
    mmseqs_pick_besthit_with_ambiguity.py \
    --m8 "${ID}.${DATABASE_NAME}.resultDB_bits_nucl_taxid.m8" \
    --ambig-frac "${AMBIG_FRAC}" \
    --out "${ID}.${DATABASE_NAME}.besthit.assigned.tsv"

    postprocess_mmseqs_taxonomy_compare.py \
    --expected-map "${ID}.mmseqs_expected.tsv" \
    --assigned "${ID}.${DATABASE_NAME}.besthit.assigned.tsv" \
    --out-summary "${ID}.${DATABASE_NAME}.bamdam_mmseqs.evaluation.summary.tsv"
    """

}

process METRICS{
    input:    
        tuple val(ID), 
        path(raw_fq), 
        path(merged_fq),
        path(rm_low_complex_fq),
        path(k2_mic_unclas_fq),
        path(mapped_bam),
        path(bamdam_bam)
        // path(filterbam),
       
    output:
        tuple val(ID), path("*.metrics")

    script:

    //set the input as optional input (if the processes are skipped)
    //def seedArg = SEED ? "--seed ${SEED}" : ""


    """
        count_fastq () {
        local f="\$1"
        [[ -s "\$f" ]] || { echo "NA"; return; }
        if gzip -t "\$f" >/dev/null 2>&1; then
            gzip -cd "\$f" 2>/dev/null | awk 'END{print (NR?NR/4:"NA")}'
        else
            awk 'END{print (NR?NR/4:"NA")}' "\$f"
        fi
        }

        count_bam () {
        local b="\$1"
        [[ -s "\$b" ]] || { echo "NA"; return; }
        if have_samtools; then
            # -F 260 = exclude unmapped (4) + secondary (256) → primary mapped alignments only
            samtools view -c -F 260 "\$b" 2>/dev/null || echo "NA"
        else
            echo "NA"
        fi

        count_raw_fq=count_fastq(raw_fq)
        count_merged_fq=count_fastq(merged_fq)
        count_rm_low_complex_fq=count_fastq(rm_low_complex_fq)
        count_k2_mic_unclas_fq=count_fastq(k2_mic_unclas_fq)
        #for different database, collect their name and read counts
        count_mapped_bam=count_bam(mapped_bam)
        count_bamdam_bam=count_bam(bamdam_bam)
        # also make the output the same format as the original metrics as input for PLOTS

        echo "$ID,\$count_raw_fq,\$count_merged_fq,\$count_rm_low_complex_fq,\$count_k2_mic_unclas_fq,\$count_mapped_bam,\$count_bamdam_bam" > $ID.metrics

    """
}

process CONCAT_METRICS {
    input: path(metrics)

    output: path("metrics")

    publishDir "${params.OUTPUT_Dir}/11_metrics", mode: "copy"

    script:

    """
    echo 'sample\tDatabase_name\traw_reads\tafter_fastp\tafter_sga_preprocess\tafter_sga_filter\tafter_prinseq\tkraken_classified\tkraken_unclassified\tmapped_phylonorway\tmapped_mito\tmapped_plastid\tmapped_mbf\tmapped_all\tafter_filterBAM\tafter_bamdam' > metrics

    cat *.metrics >> metrics


    """


}


process PLOTS{
    conda './envs/plots.yml'

    publishDir "${params.OUTPUT_Dir}/12_plots", mode: "copy"

    input:    
        tuple val(ID), 
        path(METRICS_TSV),
        path(SAMPLES_FOR_PLOTS),
        path(BAMDAM_DIR), //???----the scripts needs to change
        path(METADATA_PATH),
        path(MAP_LAST_DB_TAG),
        path(PLOT_DIR),
        val(MIN_READS),
        val(PLOTS_BAMDAM_PLOT_MODE),
        val(PLOTS_DAMAGE_THRESHOLD),
        val(PLOTS_PLOT_LOW_DAMAGE_TAXA),
        val(PLOTS_EXCLUDE_TAXA),
        val(BAMDAM_TAXA_PER_PLOT),
        val(PLOTS_LIST_TAXA_EVOLUTION_FILE),
        val(OUT_ROOT), //mmseq_dir
        val(MAP_LAST_DB_TAG),
        path(LIST_TSV),
        val(SITE_TAG),
        val(BAMDAM_MINREADS),
        val(BAMDAM_MAXDAMAGE)
            
    output:
        tuple val(ID), path("*.html"), path("*.xml")

    publishDir "${params.OUTPUT_Dir}/plots", mode: "copy"
        
    script:
    """
        100_Plots.R \
        --metrics "${METRICS_TSV}" \
        --samples "${SAMPLES_FOR_PLOTS}}" \
        --bamdam_dir "${BAMDAM_DIR}" \
        --metadata "${METADATA_PATH}" \
        --db_tag "${MAP_LAST_DB_TAG}" \
        --outdir "${PLOT_DIR}" \
        --min_reads "${MIN_READS}" \
        --bamdam_plot "${PLOTS_BAMDAM_PLOT_MODE}" \
        --damage_threshold "${PLOTS_DAMAGE_THRESHOLD}" \
        --plot_low_damage_taxa "${PLOTS_PLOT_LOW_DAMAGE_TAXA}" \
        --exclude_taxa "${PLOTS_EXCLUDE_TAXA}" \
        --taxa_per_plot "${BAMDAM_TAXA_PER_PLOT}" \
        --taxa_trend_file "${PLOTS_LIST_TAXA_EVOLUTION_FILE}" \
        --mmseqs_dir "${OUT_ROOT}/05_filtering/mmseq2"
        > "${MAP_LAST_DB_TAG}.R.out" 2>&1

        bamdam krona --in_tsv_list "${LIST_TSV}" --out_xml ${SITE_TAG}_all.xml --minreads "${BAMDAM_MINREADS}" --maxdamage "${BAMDAM_MAXDAMAGE}"

        ktImportXML -o ${SITE_TAG}_all.html ${SITE_TAG}_all.xml

    """
}
