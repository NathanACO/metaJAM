#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include the subworkflow
include { FASTP; SGA; PRINSEQ; KRAKEN2; BOWTIE2; MERGE_BAM;  MASK_REGIONS; FILTERBAM; NGSLCA;  BAMDAM; MMSEQ2 ; KRONATOOLS; METRICS; CONCAT_METRICS; PLOTS; PLOTS_KRONA_BY_SITE} from './func.nf'
include { MASK_MICROBIAL_LIKE_REGION } from './MCWorkflow/main.nf'

def validate_pipeline() {
    log.info "--- Validating Pipeline Inputs ---"

    // 1. Define the execution order and their specific tool requirements
    def process_order = [
        [name: 'FASTP',    toggle: params.ENABLE_FASTP,    req: { if(!params.FASTQ_direct_path && !params.FASTQ_list_path) error "FASTP enabled but no input FASTQ provided." }],
        [name: 'SGA',      toggle: params.ENABLE_SGA,      req: { }],
        [name: 'PRINSEQ',  toggle: params.ENABLE_PRINSEQ,  req: { }],
        [name: 'KRAKEN2',  toggle: params.ENABLE_KRAKEN_GTDB, req: { if(!params.KRAKEN2_FILTER_DATABASE) error "Missing Kraken2 Database." }],
        [name: 'MAPPING',  toggle: params.ENABLE_MAPPING,  req: { if(!params.BOWTIE2_MAPPING_DBs) error "Missing Bowtie2 Mapping DBs." }],
        [name: 'NGSLCA',   toggle: params.ENABLE_NGSLCA,   req: { if((!params.KRAKEN2_FILTER_DATABASE && !params.ACC2TAXID) || !params.NAMES || !params.NODES) error "NGSLCA requires ACC2TAXID (or bowtie2 mapping database), NAMES, and NODES files." }],
        [name: 'BAMDAM',   toggle: params.ENABLE_BAMDAM,   req: { }],
        [name: 'MMSEQS2',  toggle: params.ENABLE_MMSEQS2,  req: { if(!params.MMSEQS2_DB) error "MMSEQS2 enabled but database path is missing." }],
        [name: 'METRICS',  toggle: params.ENABLE_METRICS,  req: { }],
        [name: 'PLOTS',    toggle: params.ENABLE_PLOTS,    req: { if(!params.metadata) error "Plots enabled but 'metadata' file is missing." }]
    ]

    // 2. Identify the Start Point (The first "enable" in the list)
    def first_idx = process_order.findIndexOf { it.toggle == "enable" }
    
    if (first_idx == -1) {
        error "No modules enabled and no override files for that process was found. Nothing to do."
    }

    def startPoint = process_order[first_idx].name
    log.info "Workflow Entry Point: $startPoint"

    // 3. START POINT VALIDATION: Ensure we have the necessary input files to "jump in"
    // We only check the override for the tool that is actually starting the chain
    switch(startPoint) {
        case 'KRAKEN2': if(!params.OVERRIDE_PREPROCESSED) error "Starting at Kraken2: Provide 'OVERRIDE_PREPROCESSED'."; break
        case 'MAPPING': if(!params.OVERRIDE_LIST_KRAKEN)  error "Starting at Mapping: Provide 'OVERRIDE_LIST_KRAKEN'."; break
        case 'NGSLCA':  if(!params.OVERRIDE_LIST_BAM)     error "Starting at NGSLCA: Provide 'OVERRIDE_LIST_BAM'."; break
        case 'BAMDAM':  if(!params.OVERRIDE_LIST_NGSLCA)  error "Starting at BAMDAM: Provide 'OVERRIDE_LIST_NGSLCA'."; break
        case 'MMSEQS2': if(!params.OVERRIDE_LIST_BAMDAM)  error "Starting at MMSEQS2: Provide 'OVERRIDE_LIST_BAMDAM'."; break
        case 'METRICS': if(!params.OVERRIDE_LIST_METRICS) error "Starting at Metrics: Provide 'OVERRIDE_LIST_METRICS'."; break
    }

    // 4. SEQUENTIAL VALIDATION: Check requirements for all tools from startPoint to end
    process_order.eachWithIndex { step, i ->
        if (i >= first_idx && step.toggle == "enable") {
            step.req() 
        }
    }

    log.info "--- Validation Successful ---"
    return startPoint
}

// Global variable for the workflow block to use
workflow_start_node = validate_pipeline()


workflow {
	// input =================

	ch_sample_ids = Channel
    .fromPath(params.metadata)
    .splitCsv(header: true, sep: '\t', strip: true)
    .map { row -> row.sample }
    .unique()
	//ch_sample_ids.view()

    paired_reads = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE1", params.METAJAM_DIR+"/assets/NO_FILE2"] }
	//paired_reads = ch_sample_ids.map { id -> [id, [], []] }
	fastp_ch = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE3"] }
    preprocessed_reads = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE4"] } // for nextflow evaluation
	
	// Channel preprocessed_reads
	if (params.ENABLE_PREPROCESS == "enable" ) {
		// if (params.ENABLE_FASTP == "enable") { FASTP( input ) }

		if ( params.ENABLE_FASTP == "enable" ) {
					// allows to read fastq from a list of path
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
									.replaceAll(params.suffix_OVERRIDE_LIST, '')
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
					.map { id, reads -> 
					tuple(id.replaceAll(params.suffix_OVERRIDE_LIST, ''), reads[0], reads[1]) 
					}
			} else {
				paired_reads2 = Channel.empty()
			}

			// use both ways of specified fastq and remove duplicate
			paired_reads1.concat(paired_reads2).unique()
			.set{paired_reads}

			input_fastq = paired_reads.map { tuple(it[0], it[1], it[2], params.FASTP_OVERLAP_LEN_REQUIRE, params.FASTP_MIN_LENGTH) }
			FASTP(input_fastq)
			FASTP.out.set{fastp_ch}
		} else {
			Channel.fromPath(params.OVERRIDE_LIST_FASTP)
			.splitText()
			.map { it.trim() }
			.filter { it }
			.map { f ->
				def fq   = file(f)
				def id   = fq.baseName
					.replaceAll(/(_R?[12](_001)?|_[12])\.f(ast)?q(\.gz)?$/, '')
					.replaceAll(params.suffix_OVERRIDE_LIST, '')
				tuple(id, fq) }
			.set{fastp_ch}
		}

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
	} else if (params.ENABLE_PREPROCESS=="disable") {
		fastp_ch = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE5"] }

		preprocessed_reads = Channel
			.fromPath(params.OVERRIDE_PREPROCESSED)
			.splitText()
			.map { it.trim() }
			.filter { it }
			.map { f ->
				def fq   = file(f)
				def id   = fq.baseName
					.replaceAll(/(_R?[12](_001)?|_[12])\.f(ast)?q(\.gz)?$/, '')
					.replaceAll(params.suffix_OVERRIDE_LIST, '')
				tuple(id, fq)
			}
	}
	// preprocessed_reads.view()

	kraken_out = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE6"] }
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

	bowtie2_out = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE7"] }
	if (params.ENABLE_MAPPING == "enable") {
			Channel.fromPath( params.BOWTIE2_MAPPING_DBs )
			.splitText { it.strip( ) }
			.map { it -> 
			def name = it.tokenize('/')[-1]   // get basename
			tuple(name, file("${it}*bt2*"))}
			.groupTuple()
			.map { idx, idxs -> tuple(idx, idxs[0]) }
			.set { mapping_indexes }
			// mapping_indexes.view()

			kraken_out
			.map { it -> tuple(it[0], it[1], params.BOWTIE2_N_ALLOW_MULTIMAPPER) }
			.combine( mapping_indexes )
			.set{input_bowtie2}
			// input_bowtie2.view()
			BOWTIE2( input_bowtie2 )

			BOWTIE2.out.set { bowtie2_out }
	} else if (params.ENABLE_MAPPING=="disable") {

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
	} else if (params.ENABLE_NGSLCA != "enable" && params.OVERRIDE_LIST_NGSLCA != "") {
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
			
	} else if (params.ENABLE_BAMDAM != "enable" && params.OVERRIDE_LIST_BAMDAM != "") {
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
		MMSEQ2.out.evaluation.set{mmseq2_evaluation}
	} else {
		mmseq2_evaluation = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE8"] }
	}
	//mmseq2_evaluation.view()

	if (params.ENABLE_KRONATOOLS == "enable") { 
		if (params.ENABLE_BAMDAM != "enable" && params.OVERRIDE_LIST_BAMDAM_XML != "") {
			
			Channel.fromPath(params.OVERRIDE_LIST_BAMDAM_XML)
			.splitCsv(header: false, sep: '\t', strip: true)
			.map { row -> tuple( row[0], file(row[1]))}
			.set{bamdam_xml}
		}
		KRONATOOLS( bamdam_xml )
	}

	if (params.ENABLE_METRICS == "enable") {

		// paired_reads.view { "DEBUG paired_reads: $it" }
		// fastp_ch.view { "DEBUG fastp_ch: $it" }
		// preprocessed_reads.view { "DEBUG preprocessed_reads: $it" }
		// kraken_out.view { "DEBUG kraken_out: $it" }
		// mapped_bam.view { "DEBUG mapped_bam: $it" }
		// bamdam_bam_lca.map{it -> tuple(it[0], it[1])}.view { "DEBUG bamdam ext: $it" }
		// bamdam_bam_lca.view { "DEBUG bamdam: $it" }

		paired_reads
		.combine( fastp_ch ,by:0 )
		.combine( preprocessed_reads ,by:0 )
		.combine( kraken_out, by:0 )
		.combine( mapped_bam, by:0 )
		.combine( bamdam_bam_lca.map{it -> tuple(it[0], it[1])}, by:0 )
		.set{ input_metrics }
		input_metrics.view{ "DEBUG input_metrics: $it" }

		METRICS( input_metrics )
		// METRICS.out.collect().view()
		CONCAT_METRICS( METRICS.out.collect() )
	}

	if (params.ENABLE_PLOTS == "enable") {
		CONCAT_METRICS.out
		.combine(bamdam_bam_lca.map(it -> it[2]).collect())
		.combine(mmseq2_evaluation.collect())
		.map(it -> tuple(it[0], it[1], it[2],
			params.metadata,
			params.PLOTS_MIN_READS,
			params.PLOTS_MODE,
			params.PLOTS_DAMAGE_THRESHOLD,
			params.PLOTS_PLOT_LOW_DAMAGE_TAXA,
			params.PLOTS_EXCLUDE_TAXA,
			params.PLOTS_TAXA_PER_PLOT,
			params.PLOTS_LIST_TAXA_EVOLUTION_FILE,
			params.MAP_LAST_DB_TAG,
			params.BAMDAM_MINREADS,
			params.BAMDAM_MAXDAMAGE
		)) // params.SAMPLES_FOR_PLOTS,		
		.set{ input_plots }

		PLOTS( input_plots )

		PLOTS.out.view()

		PLOTS_KRONA_BY_SITE(
			bamdam_bam_lca.map(it -> it[3]).collect(),
			params.metadata,
			params.BAMDAM_MINREADS,
			params.BAMDAM_MAXDAMAGE
		)

	}
}