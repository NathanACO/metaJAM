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

    errorStrategy 'retry' // duplicate removal could oom without nextflow registering it as oom

    conda 'bioconda::prinseq'
    input:
    tuple val(ID), path(merged_reads), val(PRINSEQ_COMPLEXITY_METHOD), val(PRINSEQ_COMPLEXITY_THRESHOLD), val(PRINSEQ_MIN_LEN), val(PRINSEQ_DEREP)

    output:
    tuple val(ID), path("*_merged.complexity_filtered.duplicatesremoved.fastq.gz"), emit: "rm_low_complexity"    
    tuple path("*_merged.complexity_filtered.fastq.gz"), path("*_merged.low_complexity.fastq.gz")

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

        echo "low complexity filter finished"

		# 3) Remove exact fwd/rev duplicates
		prinseq-lite.pl \
		-fastq "\${GOOD}.fastq" \
		-out_good "\${DEREP_OUT}" \
		-out_bad null \
		-min_len "${PRINSEQ_MIN_LEN}" \
		-derep "${PRINSEQ_DEREP}" \
		-line_width 0

        echo "duplicate removal finished"

		gzip -f *.fastq
		"""
}

process KRAKEN2 {
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
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

output:
    tuple val(ID), path("*.bam")

    publishDir "${params.OUTPUT_Dir}/04_parallel_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper
    """
}


process GET_ACC2TAXID {
    label 'small_memory'

    conda './envs/acc2taxid.yml'

    input:
    tuple val(idx_ID), path(idxs)

    output:
    path "*.acc2taxid"

    publishDir "${params.OUTPUT_Dir}/acc2taxid", mode: "copy"

    script:
    """
        bowtie2-inspect --names $idx_ID | cut -f 1,2,3 -d" " > ${idx_ID}.contigs

        cut -f 2,3 -d' ' ${idx_ID}.contigs | sort | uniq > species

        get_acc2taxid.py species > species_taxid

        add_taxid.py ${idx_ID}.contigs species_taxid > ${idx_ID}.acc2taxid

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

        samtools sort -n -@ "${task.cpus}" -o ${ID}_merged.sorted.bam   ${ID}.merged.bam

        samtools quickcheck "${ID}_merged.sorted.bam" || { echo "Failed at merging bams of $ID"; exit 1; }
    """
}

process CONCATENATE_BEDFILES {
    label 'little_memory'

    input: path(bed_files)

    output:
        path("concatenated.bed")

    publishDir "${params.OUTPUT_Dir}/06_masked_bam", mode: "copy"
        
    script:
    """
        cat *.bed | sort -k1,1 -k2,2n > concatenated.bed
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
        samtools quickcheck "\$(basename "$bam" .bam).mic_masked.bam" || { echo "Failed at masking $bam"; exit 1; }
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
    errorStrategy = { task.exitStatus == 1 ? 'ignore' : 'retry' }
    
    input:    
        tuple val(ID), path(bam), path(NAMES), path(NODES), path(ACC2TAX)
            
    output:
        tuple val(ID), path("*.lca"), emit: lca
	path("*.log")

    publishDir "${params.OUTPUT_Dir}/07_ngslca", mode: "copy"

    errorStrategy 'ignore'
        
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

    label 'little_memory'

    publishDir "${params.OUTPUT_Dir}/08_bamdam", mode: "copy"

    input:    
        tuple val(ID), path(bam), path(lca), 
        val(STRANDED), 
        val(MINREADS), 
        val(MAXDAMAGE), 
        val(TOP_GENUS),
        val(MINCOUNT),
        val(MINSIM),
        val(MODE)
            
    output:
        tuple val(ID), path("*.small.bam"), path("*.small.lca"), path("*.tsv"), path("*.xml"), emit: lca
        tuple path("*.subs.txt"), path("*.png") // damage metircs and plots

    script:
    """
        # shrink
        bamdam shrink \
        --in_bam "${bam}" \
        --in_lca "${lca}" \
        --out_bam "${ID}.small.bam" \
        --out_lca "${ID}.small.lca" \
        --stranded "${STRANDED}" \
        --show_progress \
        --mincount ${MINCOUNT} \
        --minsim ${MINSIM}

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
        --mode ${MODE} \
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
        END{for(id in counts)printf "%d\t%s\t%s\\n",counts[id],id,name[id]}' \
        "${ID}.small.lca" \
        | sort -nrk1,1 \
        | head -n "$TOP_GENUS"
        )
  """
}

process MMSEQ2 {
    conda 'envs/mmseq2.yml'

    input:    
        tuple val(ID), path(bam), path(lca), path(lca_tsv),
        val(MIN_DMG),
        val(TOP_GENERA),
        val(MAX_READS),
        val(MIN_READS),
        val(SPACED_KMER_MODE),
        val(S),
        val(MAX_EVALUE),
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
        path("*.mmseqs_queries.fa")
        path("*.mmseqs_expected.tsv")
        path("*.besthit.assigned.tsv")
        tuple val(ID), path("*.bamdam_mmseqs.evaluation.summary.tsv"), emit: evaluation

    publishDir "${params.OUTPUT_Dir}/09_mmseq2", mode: "copy"

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
    awk -F'\t' 'NR==1 || \$12 > ${MIN_BITS}' "${DATABASE_NAME}.resultDB.m8" > "${DATABASE_NAME}.resultDB_bits.m8"

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

process METRICS {
    label 'little_memory'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), 
            path(raw_fq1), 
            path(raw_fq2), 
            path(merged_fq),
            path(rm_low_complex_fq),
            path(k2_mic_unclas_fq),
            path(mapped_bam),
            path(bamdam_bam)

        // path(filterbam),

    output:
        path("*.metrics")

    script:

    //set the input as optional input (if the processes are skipped)
    //def seedArg = SEED ? "--seed ${SEED}" : ""


    """
        count_fastq () {
        local f="\$1"

        if [ "\$f" == "NO_FILE" ]; then
            echo "NA"
            return 0
        fi

        [[ -s "\$f" ]] || { echo "NA"; return; }
        if gzip -t "\$f" >/dev/null 2>&1; then
            gzip -cd "\$f" 2>/dev/null | awk 'END{print (NR?NR/4:"NA")}'
        else
            awk 'END{print (NR?NR/4:"NA")}' "\$f"
        fi
        }

        count_bam () {
        local b="\$1"

        if [ "\$b" == "NO_FILE" ]; then
            echo "NA"
            return 0
        fi

        [[ -s "\$b" ]] || { echo "NA"; return; }
        # -F 260 = exclude unmapped (4) + secondary (256) → primary mapped alignments only
        samtools view -c -F 260 "\$b"
        }

        count_raw_fq1=\$(count_fastq $raw_fq1)
        count_raw_fq2=\$(count_fastq $raw_fq2)
        count_merged_fq=\$(count_fastq $merged_fq)
        count_rm_low_complex_fq=\$(count_fastq $rm_low_complex_fq)
        count_k2_mic_unclas_fq=\$(count_fastq $k2_mic_unclas_fq)
        #for different database, collect their name and read counts
        bam_header=""
        bam_values=""

	for bam in \$(ls *.bam | grep -v small); do
		count=\$(count_bam "\$bam")

		# build header
		bam_header=\$(echo "\${bam_header}\t\${bam}" | sed "s/${ID}.//")
    
             	# build values
                bam_values="\${bam_values}\t\${count}"
    done

		
    # build header
    count_bamdam=\$(count_bam "$bamdam_bam")
    
    # Write header
    echo -e "sample\tcount_raw_fq1\tcount_raw_fq2\tcount_merged_fq\tcount_rm_low_complex_fq\tcount_k2_mic_unclas_fq\${bam_header}\tbamdam_bam" > ${ID}.metrics
    
    # Write values
    echo -e "${ID}\t\${count_raw_fq1}\t\${count_raw_fq2}\t\${count_merged_fq}\t\${count_rm_low_complex_fq}\t\${count_k2_mic_unclas_fq}\${bam_values}\t\${count_bamdam}" >> ${ID}.metrics
    """
}

process CONCAT_METRICS {

    cpus 2
    memory 1.6.GB
    time 20.m

    input: path(metrics)

    output: path("metrics")

    publishDir "${params.OUTPUT_Dir}/10_metrics", mode: "copy"

    script:

    """
	awk 'NF && !seen[\$0]++' *.metrics > metrics

    """


}


process PLOTS {
    conda './envs/plots.yml'

    label 'little_memory'

    publishDir "${params.OUTPUT_Dir}/11_plots", mode: "copy"

    input:    
        tuple path(METRICS_TSV),
        path(BAMDAM_LCA), 
        path(mmseq2_evaluation),
        path(METADATA_PATH),
	    path(SAMPLES_FOR_PLOTS),
        val(MIN_READS),
        val(PLOTS_MODE),
        val(PLOTS_DAMAGE_THRESHOLD),
        val(PLOTS_PLOT_LOW_DAMAGE_TAXA),
        val(PLOTS_EXCLUDE_TAXA),
        val(PLOTS_TAXA_PER_PLOT),
        val(PLOTS_LIST_TAXA_EVOLUTION_FILE),
        val(MAP_LAST_DB_TAG)
            
    output:
        path("*.pdf")

    publishDir "${params.OUTPUT_Dir}/plots", mode: "copy"
        
    script:
    """
    #plots all sampoles in SAMPLES_FOR_PLOTS unless it's empty (plot all samples in metadata)
    if [[ -s "$SAMPLES_FOR_PLOTS" ]]; then
        all_samples="$SAMPLES_FOR_PLOTS"
    else
        all_samples="$METADATA_PATH"
    fi

    # collect unique sites
    sites=\$(awk -F'\t' '
    NR==1 {
        for(i=1;i<=NF;i++) if(\$i=="site") c=i;
        next
    }
    !a[\$c]++ {print \$c}
    ' "$METADATA_PATH")

    # create per-site metadata files (with header + full rows)
    for site in \$sites; do
        head -1 "\$all_samples" > "\${site}.metadata.txt"
        awk -F'\t' -v site="\$site" '
        NR==1 {
            for (i=1; i<=NF; i++) if (\$i=="site") col=i
            next
        }
        \$col==site
        ' "\$all_samples" >> "\${site}.metadata.txt"
    done

    # run plotting for each sample list / metadata file
    for samples in \$all_samples *.metadata.txt; do
        100_Plots.R \
            --metrics "$METRICS_TSV" \
            --samples "\$samples" \
            --bamdam_dir "." \
            --metadata "$METADATA_PATH" \
            --db_tag "$MAP_LAST_DB_TAG" \
            --outdir "." \
            --min_reads "$MIN_READS" \
            --bamdam_plot "$PLOTS_MODE" \
            --damage_threshold "$PLOTS_DAMAGE_THRESHOLD" \
            --plot_low_damage_taxa "$PLOTS_PLOT_LOW_DAMAGE_TAXA" \
            --exclude_taxa "$PLOTS_EXCLUDE_TAXA" \
            --taxa_per_plot "$PLOTS_TAXA_PER_PLOT" \
            --taxa_trend_file "$PLOTS_LIST_TAXA_EVOLUTION_FILE" \
            --mmseqs_dir "." \
            > "${MAP_LAST_DB_TAG}.R.out" #2>&1
    done
    """
}

process PLOTS_KRONA_BY_SITE {
    conda './envs/bamdam_kronatools.yml'

    label 'little_memory'

    publishDir "${params.OUTPUT_Dir}/11_plots", mode: "copy"

    errorStrategy { task.exitStatus == 255 ? 'ignore' : 'retry' } 

    input:
        path(tsv)
        path(METADATA_PATH)
        val(BAMDAM_MINREADS)
        val(BAMDAM_MAXDAMAGE)
            
    output: path("*.html")
        
    script:
    """
    #plot for each site:  
    sites=\$(awk -F'\\t' 'NR==1{for(i=1;i<=NF;i++) if(\$i=="site") c=i; next} !a[\$c]++{print \$c}' "${METADATA_PATH}")

    for site in \${sites}; do
        samples=\$(
            awk -F'\\t' -v site="\$site" '
            NR==1{
            for(i=1;i<=NF;i++){
                if(\$i=="sample") s=i
                if(\$i=="site")   c=i
            }
            if(!s || !c){
                print "ERROR: column not found" > "/dev/stderr"
                exit 1
            }
            next
            }
            \$c==site{print \$s}
            ' "${METADATA_PATH}" | sort -u
        )
        LIST_TSV="\${site}_all_tsv_list.txt"
        : > "\${LIST_TSV}"
        # Collect bamdam TSVs for each sample
        for s in \${samples}; do
            [[ -f "\${s}.tsv" ]] && echo "\${s}.tsv" >> "\${LIST_TSV}"
        done

        bamdam krona --in_tsv_list "\${LIST_TSV}" --out_xml "\${site}_all.xml" --minreads "${BAMDAM_MINREADS}" --maxdamage "${BAMDAM_MAXDAMAGE}"
        ktImportXML -o "\${site}_all.html" "\${site}_all.xml"
    done
        """
}

process PLOTS_KRONA_BY_SAMPLE {
    label 'small_memory'
    conda 'bioconda::krona'
    publishDir "${params.OUTPUT_Dir}/08_bamdam", mode: "copy"
    input:
        tuple val(ID), path(xml), path(tsv)
            
    output:
        path("*.html")

    script:
    """
        # Ensure TSV is available in working directory for XML to reference
        ktImportXML -o "${ID}.html" ${xml}
    """
}
