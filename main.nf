#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include the subworkflow
include { FASTP; SGA; PRINSEQ; KRAKEN2; BOWTIE2; MERGE_BAM;  MASK_REGIONS; FILTERBAM; NGSLCA;  BAMDAM; MMSEQ2 ; KRONATOOLS; METRICS; CONCAT_METRICS; PLOTS } from './func.nf'
include { MASK_MICROBIAL_LIKE_REGION } from './MCWorkflow/main.nf'

// --- PARAMETER VALIDATION (nf-core/eager style) ---
def check_params() {
    log.info "--- Validating Pipeline Inputs ---"

    // 1. Preprocess Logic
    if (params.ENABLE_PREPROCESS == "enable") {
        if (params.ENABLE_FASTP != "enable" && !params.OVERRIDE_LIST_FASTP) 
            error "[Config Error] Preprocess enabled & FASTP disabled: Please provide 'OVERRIDE_LIST_FASTP'."
        
        if (params.ENABLE_FASTP != "enable" && params.ENABLE_SGA != "enable" && params.ENABLE_PRINSEQ != "enable" && !params.OVERRIDE_PREPROCESSED)
            error "[Config Error] Preprocesses (fastp and prinseq/sga) are disabled: You must provide 'OVERRIDE_PREPROCESSED' to proceed."
    } else {
        if (!params.OVERRIDE_PREPROCESSED)
            error "[Config Error] Preprocess is disabled: You must provide 'OVERRIDE_PREPROCESSED' to proceed."
    }

    // 2. Kraken GTDB
    if (params.ENABLE_KRAKEN_GTDB == "enable") {
        if (!params.KRAKEN2_FILTER_DATABASE) error "[Config Error] Kraken2 enabled but 'KRAKEN2_FILTER_DATABASE' is missing."
    } else if (!params.OVERRIDE_LIST_KRAKEN) {
        error "[Config Error] Kraken2 disabled: Provide 'OVERRIDE_LIST_KRAKEN' for downstream modules."
    }

    // 3. Mapping & Masking
    if (params.ENABLE_MAPPING == "enable") {
        if (!params.BOWTIE2_MAPPING_DBs) error "[Config Error] Mapping enabled but 'BOWTIE2_MAPPING_DBs' is missing."
    } else if (!params.OVERRIDE_LIST_BAM) {
        error "[Config Error] Mapping disabled: Provide 'OVERRIDE_LIST_BAM'."
    }

    if (params.ENABLE_MASK_REGIONS == "enable") {
        if (params.ENABLE_GENERATE_BEDFILE_TO_MASK == "enable") {
            def mc_params = ['MCWORKFLOW_input_dir', 'MCWORKFLOW_input_list', 'MCWORKFLOW_pseudo_reads_file_dir']
            mc_params.each { if (!params."$it") error "[Config Error] Mask generation enabled but '$it' is missing." }
        } else if (!params.REGIONS_TO_MASK) {
            error "[Config Error] Masking enabled but no 'REGIONS_TO_MASK' file provided."
        }
    }

    // 4. NGSLCA
    if (params.ENABLE_NGSLCA == "enable") {
        if (!params.ACC2TAXID) error "[Config Error] NGSLCA enabled but 'ACC2TAXID' is missing."
        if (!params.NAMES || !params.NODES) error "[Config Error] NGSLCA requires 'NAMES' and 'NODES' files."
    } else if (!params.OVERRIDE_LIST_NGSLCA) {
        error "[Config Error] NGSLCA disabled: Provide 'OVERRIDE_LIST_NGSLCA'."
    }

    // 5. BAMDAM
    if (params.ENABLE_BAMDAM == "enable") {
        // No specific override needed here as it usually takes lca from NGSLCA
    } else if (!params.OVERRIDE_LIST_BAMDAM) {
        error "[Config Error] BAMDAM disabled: Downstream 'MMSEQS2' or 'PLOTS' require 'OVERRIDE_LIST_BAMDAM'."
    }

	// 6. MMSEQS2
    if (params.ENABLE_MMSEQS2 == "enable") {
        if (!params.MMSEQS2_DB) error "[Config Error] MMSEQS2 enabled but 'MMSEQS2_DB' path is missing."
        if (!params.MMSEQS2_GENERA_FILE) error "[Config Error] MMSEQS2 enabled but 'MMSEQS2_GENERA_FILE' is missing."
    }

    // 7. Krona & Plots
    if (params.ENABLE_KRONATOOLS == "enable") {
        if (params.ENABLE_BAMDAM != "enable" && !params.OVERRIDE_LIST_BAMDAM_XML) {
            error "[Config Error] Krona enabled & BAMDAM disabled: Provide 'OVERRIDE_LIST_BAMDAM_XML'."
        }
    }

	if (params.ENABLE_METRICS != "enable" && !params.OVERRIDE_LIST_METRICS) {
		error "[Config Error] MMSEQS2 enabled but 'MMSEQS2_DB' path is missing."
	}

    if (params.ENABLE_PLOTS == "enable") {
        if (!params.metadata) error "[Config Error] Plots enabled but 'metadata' file is missing."
    }

    log.info "--- Validation Successful. Starting Workflow ---"
}

// Call the function
check_params()

// --- WORKFLOW START ---
workflow {

	// allows to read fastq from a list of path
	if (params.ENABLE_PREPROCESS=="enable" && (params.FASTQ_list_path != "" || params.FASTQ_direct_path != "")) {
		if (params.FASTQ_list_path != ""){
			paired_reads1 = Channel
			.fromPath(params.FASTQ_list_path)
			.splitText()
			.map { it.trim() }
			.filter { it }
			.map { f ->
				def id = f.tokenize('/')[-1]
								.replaceAll(/(_R?1(_001)?|_1)\.f(ast)?q\.gz$/, '')
								.replaceAll(/(_R?2(_002)?|_1)\.f(ast)?q\.gz$/, '')
				tuple(id, file(f, checkIfExists: true))
			}
			.groupTuple()
			.map { id, reads -> tuple(id, reads[0], reads[1]) }
		} else {
			paired_reads1 = Channel.empty()
		}
		// allows to read fastq from direct path
		if (params.FASTQ_direct_path != ""){
			paired_reads2 = Channel
				.fromFilePairs(params.FASTQ_direct_path)
				.map { id, reads -> tuple(id, reads[0], reads[1]) }
		} else {
			paired_reads2 = Channel.empty()
		}

		// use both ways of specified fastq and remove duplicate
		paired_reads1.concat(paired_reads2).unique()
		.set{paired_reads}
		paired_reads.view()
	} else {
		paired_reads = Channel.empty()
	}

	if (params.ENABLE_MAPPING=="enable" && params.BOWTIE2_MAPPING_DBs != "" ) {
		Channel.fromPath( params.BOWTIE2_MAPPING_DBs )
		.splitText { it.strip( ) }
		.map { it -> 
		def name = it.tokenize('/')[-1]   // get basename
		tuple(name, file("${it}*bt2*"))}
		.groupTuple()
		.map { idx, idxs -> tuple(idx, idxs[0]) }
		.set { mapping_indexes }
		// mapping_indexes.view()
	} else {
		mapping_indexes = Channel.empty()
	}

	input_fastq = paired_reads.map { tuple(it[0], it[1], it[2], params.FASTP_OVERLAP_LEN_REQUIRE, params.FASTP_MIN_LENGTH) }

	preprocessed_reads = Channel.empty() // for nextflow evaluation
	fastp_ch = Channel.empty()
	// Channel preprocessed_reads
	if (params.ENABLE_PREPROCESS == "enable" ) {
		// if (params.ENABLE_FASTP == "enable") { FASTP( input ) }

		if ( params.ENABLE_FASTP == "enable" ) {
			FASTP(input_fastq)
			FASTP.out.set{fastp_ch}
		} else if (params.ENABLE_FASTP == "disable" && params.OVERRIDE_LIST_FASTP != "") {
			Channel.fromPath(params.OVERRIDE_LIST_FASTP)
			.splitText()
			.map { it.trim() }
			.filter { it }
			.map { f ->
				def fq   = file(f)
				def id   = fq.baseName
					.replaceAll(/(_R?[12](_001)?|_[12])\.f(ast)?q(\.gz)?$/, '')
				tuple(id, fq) }
			.set{fastp_ch}
		} else {
			    println "If skip fastp or fastp, you need to supply OVERRIDE_LIST_FASTP in config"
    			exit 1
		}

		// make filtering out low-complexity reads optional
		if (params.ENABLE_SGA == "enable" || params.ENABLE_PRINSEQ == "enable") {  
			if (params.ENABLE_SGA == "enable") {  
				SGA(fastp_ch .map { tuple(it[0], it[1], params.SGA_DUST_THRESHOLD) })
				// SGA.out.rm_low_complexity.view()
				}
			if (params.ENABLE_PRINSEQ == "enable") {  
				PRINSEQ(fastp_ch .map { tuple(it[0], it[1], params.PRINSEQ_COMPLEXITY_METHOD, params.PRINSEQ_COMPLEXITY_THRESHOLD, params.PRINSEQ_MIN_LEN, params.PRINSEQ_DEREP)}) 
				// PRINSEQ.out.rm_low_complexity.view()
				}

			if (params.ENABLE_LOW_COMPLEXITY_FILTER=="SGA"){
				SGA.out.rm_low_complexity.set { preprocessed_reads }
			} else if (params.ENABLE_LOW_COMPLEXITY_FILTER=="PRINSEQ") {
				PRINSEQ.out.rm_low_complexity.set { preprocessed_reads }
			}
		} else {
			preprocessed_reads = fastp_ch
		}
	} else if ((params.ENABLE_PREPROCESS=="disable" || (params.ENABLE_SGA == "disable" && params.ENABLE_PRINSEQ == "disable")) && params.OVERRIDE_PREPROCESSED != "") {
		// to test
		preprocessed_reads = Channel
			.fromPath(params.OVERRIDE_PREPROCESSED)
			.splitText()
			.map { it.trim() }
			.filter { it }
			.map { f ->
				def fq   = file(f)
				def id   = fq.baseName
					.replaceAll(/(_R?[12](_001)?|_[12])\.f(ast)?q(\.gz)?$/, '')
				tuple(id, fq)
			}
	}

	// preprocessed_reads.view()

	kraken_out = Channel.empty()
	if (params.ENABLE_KRAKEN_GTDB == "enable") { 
			KRAKEN2( preprocessed_reads.map { it -> tuple(it[0], it[1], params.KRAKEN2_FILTER_DATABASE) } )
			KRAKEN2.out.not_microbe
			.set { kraken_out }	
	} else if (params.ENABLE_KRAKEN_GTDB=="disable" && params.OVERRIDE_LIST_KRAKEN != "") {
	// } else {

		kraken_out = Channel
			.fromPath(params.OVERRIDE_LIST_KRAKEN)
			.splitText()
			.map { it.trim() }
			.filter { it }
			.map { f ->
				def fq   = file(f)
				def id   = fq.baseName
					.replaceAll(/(_R?[12](_001)?|_[12])\.f(ast)?q(\.gz)?$/, '')
				tuple(id, fq)
			}
	}
	// kraken_out.view()
	bowtie2_out = Channel.empty()
	if (params.ENABLE_MAPPING == "enable") { 
			kraken_out
			.map { it -> tuple(it[0], it[1], params.BOWTIE2_N_ALLOW_MULTIMAPPER) }
			.combine( mapping_indexes )
			.set{input_bowtie2}
			// input_bowtie2.view()
			BOWTIE2( input_bowtie2 )

			BOWTIE2.out.set { bowtie2_out }
	} else if (params.ENABLE_MAPPING=="disable" && params.OVERRIDE_LIST_BAM != "") {
		bowtie2_out = Channel
		.fromPath(params.OVERRIDE_LIST_BAM)
		.splitCsv(header: false, sep: '\t', strip: true)
		.map { row -> tuple( row[0], file(row[1]))}
	}

	bowtie2_out.groupTuple().set{mapped_bam}
	// mapped_bam.view()

	MERGE_BAM( mapped_bam )
	// MERGE_BAM.out.view()

	// MERGE_BAM.out.map { tuple(it[0], it[1], params.REGIONS_TO_MASK) } .view()

	if (params.ENABLE_MASK_REGIONS == "enable" ) { 
		if (params.ENABLE_GENERATE_BEDFILE_TO_MASK == "enable"){
			MASK_MICROBIAL_LIKE_REGION(
			params.MCWORKFLOW_input_dir, params.MCWORKFLOW_input_list, params.MCWORKFLOW_pseudo_reads_file_dir,
			params.MCWORKFLOW_type_of_pseudo_reads, params.MCWORKFLOW_n_allowed_multimappers,
			"${params.OUTPUT_Dir}/bedfile_for_masking", './MCWorkflow', './MCWorkflow/GTDB_fna2name.txt'
			)
			MASK_MICROBIAL_LIKE_REGION.out.set{regions_to_mask}
		} else {
			Channel.fromPath( params.REGIONS_TO_MASK ).set{regions_to_mask}
		}
		// regions_to_mask.view()

		MERGE_BAM.out
		.combine(regions_to_mask)
		.set{ input_mask }
		//input_mask.view()

		MASK_REGIONS( input_mask )
		MASK_REGIONS.out.set{ bam }

	} else {
		MERGE_BAM.out.set{ bam }
	}
	// if (params.ENABLE_FILTERBAM == "enable") { 
	// 	FILTERBAM( input )
	// 	FILTERBAM.out.set{ bam }
	// 	} else { MERGE_BAM.out.set{ bam } }

	// bam.map { tuple(it[0], it[1], params.NAMES, params.NODES, params.ACC2TAXID )}.view()

	// lca = Channel.empty()
	if (params.ENABLE_NGSLCA == "enable") { 
		
		// acc2taxid is only needed by ngsLCA
		def acc2taxid = file(params.ACC2TAXID) 
		def acc2taxid_exists = acc2taxid.exists()
		if ( !acc2taxid.exists() ) {

			Channel.fromPath( params.BOWTIE2_MAPPING_DBs )
			.splitText { it.strip( ) }
			.map { it -> 
			def name = it.tokenize('/')[-1]   // get basename
			tuple(name, file("${it}*bt2*"))}
			.groupTuple()
			.map { idx, idxs -> tuple(idx, idxs[0]) }
			.set { mapping_indexes }

			GET_ACC2TAXID.out.collectFile(name: 'acc2taxid.txt', newLine: true)
			.set{acc2taxid}
		}

		NGSLCA( bam.map { tuple(it[0], it[1], params.NAMES, params.NODES, acc2taxid )} )
		NGSLCA.out.lca.set{lca}
	} else if (params.ENABLE_NGSLCA == "disable" && params.OVERRIDE_LIST_NGSLCA != "") {
		Channel.fromPath(params.OVERRIDE_LIST_NGSLCA)
		.splitCsv(header: false, sep: '\t', strip: true)
		.map { row -> tuple( row[0], file(row[1]))}
		.set{lca}
	}
	// lca.view()

	if (params.ENABLE_BAMDAM == "enable" ) {

		MERGE_BAM.out
		.combine( lca ,by:0 )
		.map(it -> tuple(it[0], it[1], it[2],
		params.BAMDAM_STRANDED,
		params.BAMDAM_MINREADS,
		params.BAMDAM_MAXDAMAGE,
		params.BAMDAM_TOP_GENUS,
		params.BAMDAM_MINCOUNT,
        params.BAMDAM_MINSIM,
        params.BAMDAM_MODE
		))
		.set{ input_bamdam }
		//input_bamdam.view()

		BAMDAM( input_bamdam )
		BAMDAM.out.lca.set{bamdam_bam_lca}

		BAMDAM.out.xml.set{bamdam_xml}
		// bamdam_bam_lca.view()
			
	} else if (params.ENABLE_BAMDAM == "disable" && params.OVERRIDE_LIST_BAMDAM != "") {
		Channel
		.fromPath(params.OVERRIDE_LIST_BAMDAM)
		.splitCsv(header: false, sep: '\t', strip: true)
		.map { row -> 
			tuple( row[0], file(row[1]), file(row[2]), file(row[3]))}
		.set{bamdam_bam_lca}

		// bamdam_bam_lca.view()
	}

	if (params.ENABLE_MMSEQS2 == "enable") {
		MMSEQS2_GENERA_FILE = Channel.fromPath(params.MMSEQS2_GENERA_FILE, checkIfExists:true)

		bamdam_bam_lca
		.map(it -> tuple(it[0], it[1], it[2], it[3],
		params.MMSEQS2_MIN_DMG,
		params.MMSEQS2_TOP_GENERA,
		params.MMSEQS2_MAX_READS,
		params.MMSEQS2_MIN_READS,
		params.MMSEQS2_SPACED_KMER_MODE,
		params.MMSEQS2_S,
		params.MMSEQS2_MAX_EVALUE,
		params.MMSEQS2_MIN_LENGTH,
		params.MMSEQS2_MIN_SEQID,
		params.MMSEQS2_MAX_SEQS,
		params.MMSEQS2_MIN_QUERY_COV,
		params.MMSEQS2_SPLIT_MEM_LIMIT,
		params.MMSEQS2_AMBIG_FRAC,
		params.MMSEQS2_DATABASE_NAME,
		params.MMSEQS2_DB,
		params.MMSEQS2_SEED,
		params.MMSEQS2_MIN_BITS,
		params.MMSEQS2_TAXADB_SQLITE
		))
		.combine(MMSEQS2_GENERA_FILE) //optional input
		.set{ input_mmseq2 }
		//input_mmseq2.view()
		
		MMSEQ2( input_mmseq2 )
		MMSEQ2.out.evaluation.collect().set{mmseq2_evaluation}
	} else {
		mmseq2_evaluation = Channel.empty()
	}

	if (params.ENABLE_KRONATOOLS == "enable") { 
		if (params.ENABLE_BAMDAM == "disable" && params.OVERRIDE_LIST_BAMDAM_XML != "") {
			
			Channel.fromPath(params.OVERRIDE_LIST_BAMDAM_XML)
			.splitCsv(header: false, sep: '\t', strip: true)
			.map { row -> tuple( row[0], file(row[1]))}
			.set{bamdam_xml}
		}
		KRONATOOLS( bamdam_xml )
	}

	metrics = Channel.empty()
	if (params.ENABLE_METRICS == "enable") {

		paired_reads
		.combine( fastp_ch ,by:0 )
		.combine( preprocessed_reads ,by:0 )
		.combine( kraken_out, by:0 )
		.combine( mapped_bam, by:0 )
		.combine( bamdam_bam_lca.map{it -> tuple(it[0], it[1])}, by:0 )
		.set{ input_metrics }
		input_metrics.view()

		METRICS( input_metrics )
		METRICS.out.collect().view()
		CONCAT_METRICS( METRICS.out.collect() )

		CONCAT_METRICS.out.set{ metrics }

	} else {
		Channel.fromPath(params.OVERRIDE_LIST_METRICS)
		.set{ metrics }
	}

	if (params.ENABLE_PLOTS == "enable") {
		metrics
		.combine(bamdam_bam_lca.map(it -> it[2]).collect())
		.combine(mmseq2_evaluation)
		.map(it -> tuple(it[0], it[1], it[2],
			params.SAMPLES_FOR_PLOTS,		
			params.metadata,
			params.MAP_LAST_DB_TAG,
			params.PLOT_DIR,
			params.MIN_READS,
			params.PLOTS_BAMDAM_PLOT_MODE,
			params.PLOTS_DAMAGE_THRESHOLD,
			params.PLOTS_PLOT_LOW_DAMAGE_TAXA,
			params.PLOTS_EXCLUDE_TAXA,
			params.BAMDAM_TAXA_PER_PLOT,
			params.PLOTS_LIST_TAXA_EVOLUTION_FILE,
			params.MAP_LAST_DB_TAG,
			params.BAMDAM_MINREADS,
			params.BAMDAM_MAXDAMAGE
		))
		.set{ input_plots }

		PLOTS( input_plots )
	}
		
	
}