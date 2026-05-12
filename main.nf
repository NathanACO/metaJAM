nextflow.enable.dsl = 2

// include the subworkflow
include { FASTP; SGA; PRINSEQ; KRAKEN2; BOWTIE2; MERGE_BAM; CONCATENATE_BEDFILES; MASK_REGIONS; FILTERBAM; NGSLCA;  BAMDAM; MMSEQ2 ; PLOTS_KRONA_BY_SAMPLE; METRICS; CONCAT_METRICS; PLOTS; PLOTS_KRONA_BY_SITE; GET_ACC2TAXID } from './func.nf'
include { MASK_MICROBIAL_LIKE_REGION } from './subworkflows/mcworkflow.nf'
include { SEQUENTIAL_MAP } from './subworkflows/sequential_map.nf'

def validate_pipeline() {

    log.info "--- Validating Pipeline Inputs ---"

	if (!params.metadata) error "Missing required file parameter: metadata"
	if (params.METAJAM_DIR == "" || params.METAJAM_DIR == null) error "Missing required value parameter: METAJAM_DIR"

    // 1. Define the execution order and their specific tool requirements
    def process_order = [
		[name: 'FASTP',    toggle: params.ENABLE_FASTP,    req: {
			if ((params.FASTQ_direct_path == "" || params.FASTQ_direct_path == null) && !params.FASTQ_list_path) {
				error "FASTP enabled but no input FASTQ provided. Set FASTQ_direct_path or FASTQ_list_path."
			}
		}],
        [name: 'SGA',      toggle: params.ENABLE_SGA,      req: {if (params.ENABLE_SGA != "enable" && !params.OVERRIDE_PREPROCESSED) error "SGA enabled but params.OVERRIDE_PREPROCESSED is not provided."}],
        [name: 'PRINSEQ',  toggle: params.ENABLE_PRINSEQ,  req: {if (params.ENABLE_PRINSEQ != "enable" && !params.OVERRIDE_PREPROCESSED) error "PRINSEQ enabled but params.OVERRIDE_PREPROCESSED is not provided."}],
		[name: 'KRAKEN2',  toggle: params.ENABLE_KRAKEN_GTDB, req: {
			if (!params.KRAKEN2_FILTER_DATABASE) error "Missing required file parameter: KRAKEN2_FILTER_DATABASE"
		}],
		[name: 'MAPPING',  toggle: params.ENABLE_MAPPING,  req: {
			if (!params.BOWTIE2_MAPPING_DB1) error "Missing required file parameter: BOWTIE2_MAPPING_DB1"
			if (params.BOWTIE2_N_ALLOW_MULTIMAPPER == "" || params.BOWTIE2_N_ALLOW_MULTIMAPPER == null) error "Missing required value parameter: BOWTIE2_N_ALLOW_MULTIMAPPER"
			if (params.ENABLE_PARALLEL_MAPPING != "enable" && params.USE_MAPPING == "PARALLEL_MAPPING") error "params.USE_MAPPING == \"PARALLEL_MAPPING\" requires ENABLE_PARALLEL_MAPPING to be enabled."
			if (params.ENABLE_SEQUENTIAL_MAPPING != "enable" && params.USE_MAPPING == "SEQUENTIAL_MAPPING") error "params.USE_MAPPING == \"SEQUENTIAL_MAPPING\" requires ENABLE_SEQUENTIAL_MAPPING to be enabled."
			if (params.ENABLE_PARALLEL_MAPPING != "enable" && params.ENABLE_SEQUENTIAL_MAPPING != "enable") error "As mapping is enabled with params.ENABLE_MAPPING, you need to enable at least one of ENABLE_PARALLEL_MAPPING or ENABLE_SEQUENTIAL_MAPPING."
		}],
		[name: 'MASKING',  toggle: params.ENABLE_MASK_REGIONS, req: {
			if (params.ENABLE_GENERATE_BEDFILE_TO_MASK != "enable" && !params.REGIONS_TO_MASK) error "Masking enabled but required either enable running mcworkflow with ENABLE_GENERATE_BEDFILE_TO_MASK or supply REGIONS_TO_MASK file"
		}],
		[name: 'NGSLCA',   toggle: params.ENABLE_NGSLCA,   req: {
			if (!params.NAMES) error "Missing required file parameter: NAMES"
			if (!params.NODES) error "Missing required file parameter: NODES"
			if (!params.ACC2TAXID && !params.BOWTIE2_MAPPING_DB1) error "NGSLCA requires ACC2TAXID or BOWTIE2_MAPPING_DB1."
		}],
        [name: 'BAMDAM',   toggle: params.ENABLE_BAMDAM,   req: { 
			if ((params.ENABLE_MAPPING != "enable" && !params.OVERRIDE_LIST_BAM) || (params.ENABLE_NGSLCA != "enable" && !params.OVERRIDE_LIST_NGSLCA)) error "BAMDAM requires either mapping or NGSLCA to be enabled or overridden."		
		}],
		[name: 'MMSEQS2',  toggle: params.ENABLE_MMSEQS2,  req: {
			if (!params.MMSEQS2_DB) error "Missing required file parameter: MMSEQS2_DB"
			if (!params.MMSEQS2_TAXADB_SQLITE) error "Missing required file parameter: MMSEQS2_TAXADB_SQLITE"
			if (!params.MMSEQS2_GENERA_FILE) error "Missing required file parameter: MMSEQS2_GENERA_FILE"
			if (params.MMSEQS2_DATABASE_NAME == "" || params.MMSEQS2_DATABASE_NAME == null) error "Missing required value parameter: MMSEQS2_DATABASE_NAME"
		}],
        [name: 'METRICS',  toggle: params.ENABLE_METRICS,  req: { }],
		[name: 'PLOTS',    toggle: params.ENABLE_PLOTS,    req: {
			if (!params.metadata) error "Plots enabled but required file parameter metadata is missing."
			if (params.MAP_LAST_DB_TAG == "" || params.MAP_LAST_DB_TAG == null) error "Missing required value parameter: MAP_LAST_DB_TAG"
		}],
		[name: 'PLOTS_KRONA_BY_SAMPLE', toggle: params.ENABLE_KRONATOOLS, req: { if(params.ENABLE_BAMDAM != "enable" && !params.OVERRIDE_LIST_BAMDAM) error "Kronatools enabled without BamDam requires OVERRIDE_LIST_BAMDAM file parameter." }],
		[name: 'PLOTS_KRONA_BY_SITE', toggle: params.ENABLE_PLOTS, req: {
			if (params.ENABLE_BAMDAM != "enable" && !params.OVERRIDE_LIST_BAMDAM) error "Plots:PLOTS_KRONA_BY_SITE enabled but neither upstream bamdam is enabled nor params.OVERRIDE_LIST_BAMDAM is provided. Please change of the settings."
			if (!params.metadata) error "Krona plots enabled but required file parameter metadata is missing."
			if (params.BAMDAM_MINREADS == "" || params.BAMDAM_MINREADS == null) error "Missing required value parameter: BAMDAM_MINREADS"
			if (params.BAMDAM_MAXDAMAGE == "" || params.BAMDAM_MAXDAMAGE == null) error "Missing required value parameter: BAMDAM_MAXDAMAGE"
			}
    	],
		[name: 'MCworkflow',  toggle: params.ENABLE_GENERATE_BEDFILE_TO_MASK, req: {
		if (!params.BOWTIE2_MAPPING_DB1 && (!params.MCWORKFLOW_input_dir && !params.MCWORKFLOW_input_list)) error "Missing required file parameter: MCWORKFLOW_input_list/MCWORKFLOW_input_dir or BOWTIE2_MAPPING_DB1. At least one of them is required to generate mapping indexes for masking."
		if (!params.MCWORKFLOW_pseudo_reads_file_dir) error "Missing required file parameter: MCWORKFLOW_pseudo_reads_file_dir"
		if (params.MCWORKFLOW_type_of_pseudo_reads == "" || params.MCWORKFLOW_type_of_pseudo_reads == null) error "Missing required value parameter: MCWORKFLOW_type_of_pseudo_reads"
		if (params.MCWORKFLOW_n_allowed_multimappers == "" || params.MCWORKFLOW_n_allowed_multimappers == null) error "Missing required value parameter: MCWORKFLOW_n_allowed_multimappers"

		}]
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
		case 'MASKING': if(!params.OVERRIDE_LIST_BAM)  error "Starting at Masking: Provide 'OVERRIDE_LIST_BAM'."; break
        case 'NGSLCA':  if(!params.OVERRIDE_LIST_BAM)     error "Starting at NGSLCA: Provide 'OVERRIDE_LIST_BAM'."; break
        case 'BAMDAM':  if(!params.OVERRIDE_LIST_NGSLCA)  error "Starting at BAMDAM: Provide 'OVERRIDE_LIST_NGSLCA'."; break
        case 'MMSEQS2': if(!params.OVERRIDE_LIST_BAMDAM)  error "Starting at MMSEQS2: Provide 'OVERRIDE_LIST_BAMDAM'."; break
		case 'PLOTS': if(!params.OVERRIDE_LIST_METRICS) error "Starting at Plots: Provide 'OVERRIDE_LIST_METRICS'."; break
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

workflow_entry_point = validate_pipeline()

// --- WORKFLOW START ---
workflow {

		// use MCworkflow to generate bedfile for microbial-like regions
	if (params.ENABLE_GENERATE_BEDFILE_TO_MASK == "enable"){
		MASK_MICROBIAL_LIKE_REGION( params.MCWORKFLOW_input_dir, params.MCWORKFLOW_input_list,
		params.BOWTIE2_MAPPING_DB1,params.BOWTIE2_MAPPING_DB2,params.BOWTIE2_MAPPING_DB3,
		params.BOWTIE2_MAPPING_DB4,params.BOWTIE2_MAPPING_DB5,params.BOWTIE2_MAPPING_DB6,
		params.BOWTIE2_MAPPING_DB7,params.BOWTIE2_MAPPING_DB8,params.BOWTIE2_MAPPING_DB9,
		params.BOWTIE2_MAPPING_DB10, 
		params.MCWORKFLOW_pseudo_reads_file_dir,
		params.MCWORKFLOW_type_of_pseudo_reads, params.MCWORKFLOW_n_allowed_multimappers,
		"${params.METAJAM_DIR}", "${params.METAJAM_DIR}/assets/GTDB_fna2name.txt", params.ENABLE_MASK_FASTA
		)

		CONCATENATE_BEDFILES( MASK_MICROBIAL_LIKE_REGION.out.map{ it -> it[1]}.collect() )

		CONCATENATE_BEDFILES.out.set{regions_to_mask}

	} else {
		Channel.fromPath( params.REGIONS_TO_MASK ).set{regions_to_mask}
	}

	// run the main workflow only if at least one of the main steps is enabled
	if (
	params.ENABLE_PREPROCESS == "enable" || params.ENABLE_KRAKEN_GTDB == "enable" ||
    params.ENABLE_MAPPING == "enable" || params.ENABLE_MASK_REGIONS == "enable" ||
    params.ENABLE_NGSLCA == "enable" || params.ENABLE_BAMDAM == "enable" ||
    params.ENABLE_KRONATOOLS == "enable" || params.ENABLE_MMSEQS2 == "enable" ||
    params.ENABLE_METRICS == "enable" || params.ENABLE_PLOTS == "enable"
	) 
	{

		ch_sample_ids = Channel
		.fromPath(params.metadata)
		.splitCsv(header: true, sep: '\t', strip: true)
		.map { row -> row.sample }
		.unique()
		//ch_sample_ids.view()

		// Channel preprocessed_reads
		if (params.ENABLE_PREPROCESS == "enable" ) {
			// if (params.ENABLE_FASTP == "enable") { FASTP( input ) }

			if (params.FASTQ_list_path || params.FASTQ_direct_path) {
				
				// allows to read fastq from a list of path
				if (params.FASTQ_list_path){
					paired_reads1 = Channel
					.fromPath(params.FASTQ_list_path)
					.splitCsv(header: false, sep: '\t', strip: true)
					.map { row -> tuple( row[0], file(row[1]), file(row[2]))}
					// .groupTuple()
					// .map { id, reads -> tuple(id, reads[0], reads[1]) }
					//paired_reads1.view{"Debug paired_reads1 from FASTQ_list_path: ${it}"}

				} else {
					paired_reads1 = Channel.empty()
				}
				// allows to read fastq from direct path
				if (params.FASTQ_direct_path != ""){
					paired_reads2 = Channel
						.fromFilePairs(params.FASTQ_direct_path)
						.map { id, reads -> tuple(id.replaceAll(params.suffix_OVERRIDE_LIST, ''), reads[0], reads[1]) }
				} else {
					paired_reads2 = Channel.empty()
				}

				// use both ways of specified fastq and remove duplicate
				paired_reads1.concat(paired_reads2).unique()
				.set{paired_reads}
			} else {
				paired_reads = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE1", params.METAJAM_DIR+"/assets/NO_FILE2"] }
			}
			
			//paired_reads.view()

			if ( params.ENABLE_FASTP == "enable" ) {
				input_fastq = paired_reads.map { tuple(it[0], it[1], it[2], params.FASTP_OVERLAP_LEN_REQUIRE, params.FASTP_MIN_LENGTH) }
				FASTP(input_fastq)
				FASTP.out.set{fastp_ch1}
			} else {fastp_ch1 = Channel.empty()}
			
			if (workflow_entry_point == "SGA" || workflow_entry_point == "PRINSEQ" || params.OVERRIDE_LIST_FASTP ) {
				Channel
				.fromPath(params.OVERRIDE_LIST_FASTP)
				.splitCsv(header: false, sep: '\t', strip: true)
				.map { row -> tuple( row[0], file(row[1]))}
				.set{fastp_ch2}
			} else {fastp_ch2 = Channel.empty()}

			fastp_ch1.concat(fastp_ch2).unique().set{fastp_ch}

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
					SGA.out.rm_low_complexity.set { preprocessed_reads1 }
				} else if (params.ENABLE_LOW_COMPLEXITY_FILTER=="PRINSEQ") {
					PRINSEQ.out.rm_low_complexity.set { preprocessed_reads1 }
				}
			} else {
				preprocessed_reads1 = fastp_ch
			}
		//} else if (workflow_entry_point == "KRAKEN2" || params.OVERRIDE_PREPROCESSED) {
		} else {
			paired_reads = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE1", params.METAJAM_DIR+"/assets/NO_FILE2"] }
			fastp_ch = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE3"] }
			preprocessed_reads1 = Channel.empty()  // for nextflow evaluation
		}
		
		if (workflow_entry_point == "KRAKEN2" || params.OVERRIDE_PREPROCESSED) {
			// to test
			preprocessed_reads2 = Channel
				.fromPath(params.OVERRIDE_PREPROCESSED)
				.splitCsv(header: false, sep: '\t', strip: true)
				.map { row -> tuple( row[0], file(row[1]))}
		} else {
			preprocessed_reads2 = Channel.empty()
		}
		
		preprocessed_reads = preprocessed_reads1
		.concat(preprocessed_reads2)
		.unique()
		.ifEmpty {
			ch_sample_ids.map { id -> [id, params.METAJAM_DIR + "/assets/NO_FILE4"] }
		}

		// preprocessed_reads.view()
		if (params.ENABLE_KRAKEN_GTDB == "enable" || workflow_entry_point == "MAPPING" || params.OVERRIDE_LIST_KRAKEN ) {
			if (params.ENABLE_KRAKEN_GTDB == "enable") { 
					KRAKEN2( preprocessed_reads.map { it -> tuple(it[0], it[1], params.KRAKEN2_FILTER_DATABASE) } )
					KRAKEN2.out.not_microbe
					.set { kraken_out1 }	
			} else {
				kraken_out1 = Channel.empty()
			} 
			
			if (workflow_entry_point == "MAPPING" || params.OVERRIDE_LIST_KRAKEN ) {

				kraken_out2 = Channel
					.fromPath(params.OVERRIDE_LIST_KRAKEN)
					.splitCsv(header: false, sep: '\t', strip: true)
					.map { row -> tuple( row[0], file(row[1]))}
			} else {
				kraken_out2 = Channel.empty()
			}
			//kraken_out2.view()
			kraken_out1.concat(kraken_out2).unique().set{kraken_out}

			kraken_out.set{kraken_out_metrics}

		} else {
			preprocessed_reads.set{kraken_out}
			kraken_out_metrics = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE5"] }
		}
		//kraken_out.view()

		bowtie2_out1 = Channel.empty()
		if (params.ENABLE_MAPPING == "enable") { 

				// Channel.fromPath( params.BOWTIE2_MAPPING_DBs )
				// .splitText { it.strip( ) }
				// .map { it -> 
				// def name = it.tokenize('/')[-1]   // get basename
				// tuple(name, file("${it}*.bt2*"))}
				// .groupTuple()
				// .map { idx, idxs -> tuple(idx, idxs[0]) }
				// .set { mapping_indexes }

				// read mapping indexes for debugging

				if (params.ENABLE_MAPPING == "enable") {

					if (params.ENABLE_PARALLEL_MAPPING == "enable") {

						mapping_index1  = params.BOWTIE2_MAPPING_DB1  ? Channel.value(params.BOWTIE2_MAPPING_DB1).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index2  = params.BOWTIE2_MAPPING_DB2  ? Channel.value(params.BOWTIE2_MAPPING_DB2).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index3  = params.BOWTIE2_MAPPING_DB3  ? Channel.value(params.BOWTIE2_MAPPING_DB3).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index4  = params.BOWTIE2_MAPPING_DB4  ? Channel.value(params.BOWTIE2_MAPPING_DB4).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index5  = params.BOWTIE2_MAPPING_DB5  ? Channel.value(params.BOWTIE2_MAPPING_DB5).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index6  = params.BOWTIE2_MAPPING_DB6  ? Channel.value(params.BOWTIE2_MAPPING_DB6).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index7  = params.BOWTIE2_MAPPING_DB7  ? Channel.value(params.BOWTIE2_MAPPING_DB7).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index8  = params.BOWTIE2_MAPPING_DB8  ? Channel.value(params.BOWTIE2_MAPPING_DB8).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index9  = params.BOWTIE2_MAPPING_DB9  ? Channel.value(params.BOWTIE2_MAPPING_DB9).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }  : Channel.empty()
						mapping_index10 = params.BOWTIE2_MAPPING_DB10 ? Channel.value(params.BOWTIE2_MAPPING_DB10).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) } : Channel.empty()

						mapping_indexes= Channel.empty()
						.concat(mapping_index1).concat(mapping_index2).concat(mapping_index3)
						.concat(mapping_index4).concat(mapping_index5).concat(mapping_index6)
						.concat(mapping_index7).concat(mapping_index8).concat(mapping_index9)
						.concat(mapping_index10)

						// mapping_indexes.view{"Debug mapping_indexes in parallel mapping: ${it}"}

						kraken_out
						.map { it -> tuple(it[0], it[1], params.BOWTIE2_N_ALLOW_MULTIMAPPER) }
						.combine( mapping_indexes )
						.set{input_bowtie2}
						// input_bowtie2.view()
						BOWTIE2( input_bowtie2 )
						if (params.USE_MAPPING == "PARALLEL_MAPPING") { BOWTIE2.out.set { bowtie2_out1 } }
					} 

					if (params.ENABLE_SEQUENTIAL_MAPPING == "enable") {
						SEQUENTIAL_MAP( kraken_out, params.BOWTIE2_N_ALLOW_MULTIMAPPER )
						
						if (params.USE_MAPPING == "SEQUENTIAL_MAPPING") { SEQUENTIAL_MAP.out.set { bowtie2_out1 } }
					} 
				}
		} //else {bowtie2_out1 = Channel.empty()} 
		
		bowtie2_out2 = Channel.empty()
		if (workflow_entry_point == "NGSLCA" || workflow_entry_point == "MASKING" || params.OVERRIDE_LIST_BAM ) {
			if (workflow_entry_point == "NGSLCA"){mapping_indexes = Channel.empty()}
			bowtie2_out2 = Channel
			.fromPath(params.OVERRIDE_LIST_BAM)
			.splitCsv(header: false, sep: '\t', strip: true)
			.map { row -> tuple( row[0], file(row[1]))}
		} //else {bowtie2_out2 = Channel.empty()}

		// bowtie2_out1.concat(bowtie2_out2).unique().set{bowtie2_out}

		bowtie2_out = bowtie2_out1.concat(bowtie2_out2).unique()
		.map{id, bams -> tuple(id, bams)}
		.ifEmpty {
			ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE6"] }
		}
		// bowtie2_out.view()

		bowtie2_out
		.map { id, files -> tuple(id, files) }
		.groupTuple()
		.set{mapped_bam}
		mapped_bam.view()


		MERGE_BAM( mapped_bam )
		MERGE_BAM.out.set{ merged_bam }

		// regions_to_mask.view()

		// microbial-like regions
		if (params.ENABLE_MASK_REGIONS == "enable" ) { 

			merged_bam
			.combine(regions_to_mask)
			.set{ input_mask }
			//input_mask.view()

			MASK_REGIONS( input_mask )
			MASK_REGIONS.out.set{ bam }

		} else {
			merged_bam.set{ bam }
		}
		// if (params.ENABLE_FILTERBAM == "enable") { 
		// 	FILTERBAM( input )
		// 	FILTERBAM.out.set{ bam }
		// 	} else { merged_bam.set{ bam } }

		// bam.map { tuple(it[0], it[1], params.NAMES, params.NODES, params.ACC2TAXID )}.view()

		// lca = Channel.empty()
		lca1 = Channel.empty()
		if (params.ENABLE_NGSLCA == "enable") { 
			
			// acc2taxid is only needed by ngsLCA
			def acc2taxid = file(params.ACC2TAXID) 
			def acc2taxid_exists = acc2taxid.exists()
			if ( !acc2taxid.exists() ) {
				mapping_indexes | GET_ACC2TAXID | collectFile // it requires species name just after contig name in ref, also could generate NA for some species
			}

			NGSLCA( bam.map { tuple(it[0], it[1], params.NAMES, params.NODES, acc2taxid )} )
			NGSLCA.out.lca.set{lca1}
		}
		
		lca2 = Channel.empty()
		if (workflow_entry_point == "BAMDAM" || params.OVERRIDE_LIST_NGSLCA ) {
			Channel.fromPath(params.OVERRIDE_LIST_NGSLCA)
			.splitCsv(header: false, sep: '\t', strip: true)
			.map { row -> tuple( row[0], file(row[1]))}
			.set{lca2}
		}
		lca1.concat(lca2).unique().set{lca}
		// lca.view()

		bamdam_bam_lca2 = Channel.empty()
		bamdam_xml2 = Channel.empty()
		bamdam_bam_lca1 = Channel.empty()
		bamdam_xml1 = Channel.empty()
		if (params.ENABLE_BAMDAM == "enable" ) {

			merged_bam
			.combine( lca ,by:0 )
			.map{it -> tuple(it[0], it[1], it[2],
			params.BAMDAM_STRANDED,
			params.BAMDAM_MINREADS,
			params.BAMDAM_MAXDAMAGE,
			params.BAMDAM_TOP_GENUS,
			params.BAMDAM_MINCOUNT,
			params.BAMDAM_MINSIM,
			params.BAMDAM_MODE
			)}
			.set{ input_bamdam }
			//input_bamdam.view()

			BAMDAM( input_bamdam )
			BAMDAM.out.lca
			.map { row -> tuple( row[0], file(row[1]), file(row[2]), file(row[3]) )}
			.set{bamdam_bam_lca1}

			BAMDAM.out.lca
			.map { row -> tuple( row[0], file(row[4]), file(row[3]))}
			.set{bamdam_xml1}
			// bamdam_bam_lca.view()
				
		}
		
		if (workflow_entry_point == "MMSEQS2" || params.OVERRIDE_LIST_BAMDAM ) {
			Channel
			.fromPath(params.OVERRIDE_LIST_BAMDAM)
			.splitCsv(header: false, sep: '\t', strip: true)
			.map { row -> 
				tuple( row[0], file(row[1]), file(row[2]), file(row[3]),file(row[4]))}
			.set{bamdam_bam_override}

			bamdam_bam_override
			.map { row -> tuple( row[0], file(row[1]), file(row[2]), file(row[3]) )}
			.set{bamdam_bam_lca2}

			bamdam_bam_override
			.map { row -> tuple( row[0], file(row[4]), file(row[3]))}
			.set{bamdam_xml2}

			// bamdam_bam_lca.view()
		}

		bamdam_bam_lca1.concat(bamdam_bam_lca2).unique().set{bamdam_bam_lca}

		mmseq2_evaluation1 = Channel.empty()

		if (params.ENABLE_MMSEQS2 == "enable" || workflow_entry_point == "MMSEQS2") {

			if (params.ENABLE_MMSEQS2 == "enable") {
				MMSEQS2_GENERA_FILE = Channel.fromPath(params.MMSEQS2_GENERA_FILE, checkIfExists:true)

				bamdam_bam_lca
				.map{it -> tuple(it[0], it[1], it[2], it[3],
				params.MMSEQS2_MIN_DMG,
				params.MMSEQS2_TOP_GENERA,
				params.MMSEQS2_MAX_READS,
				params.MMSEQS2_MIN_READS,
				params.MMSEQS2_SPACED_KMER_MODE,
				params.MMSEQS2_S,
				params.MMSEQS2_MAX_EVALUE,
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
				)}
				.combine(MMSEQS2_GENERA_FILE) //optional input
				.set{ input_mmseq2 }
				//input_mmseq2.view()
				
				MMSEQ2( input_mmseq2 )
				MMSEQ2.out.evaluation.collect().set{mmseq2_evaluation1}
			}
			
			mmseq2_evaluation2 = Channel.empty()
			if ( params.OVERRIDE_LIST_MMSEQ2 ) {
				Channel
				.fromPath(params.OVERRIDE_LIST_MMSEQ2)
				.splitCsv(header: false, sep: '\t', strip: true)
				.map { row -> 
					tuple( row[0], file(row[1]))}
				.set{mmseq2_evaluation2}
			}

			mmseq2_evaluation = mmseq2_evaluation1.concat(mmseq2_evaluation2).unique()
		} else {mmseq2_evaluation = ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE7"]}}

		if (params.ENABLE_KRONATOOLS == "enable") { 

			bamdam_xml1.concat(bamdam_xml2).unique().set{bamdam_xml}

			PLOTS_KRONA_BY_SAMPLE( bamdam_xml )
		}

		metrics = Channel.empty()
		if (params.ENABLE_METRICS == "enable") {

			// added to avoid calling validations below if the workflow ran into error before or terminated
			onError:
				return
			//check if any channel does not match the ch_sample_ids in metadata so it does not hang silently
			paired_reads = paired_reads.combine( ch_sample_ids ,by:0 )
			.ifEmpty {
				log.error "paired_reads (params.FASTQ_list_path or params.FASTQ_direct_path) does not match with the sample IDs in params.metadata, or row in inputs not spaced with tabs. Please check your input files and metadata."
				paired_reads.view{log.error "Debug paired_reads: ${it}"}
				ch_sample_ids.view{log.error "Debug ch_sample_ids: ${it}"}
				System.exit(1)
			}

			fastp_ch = fastp_ch.combine(ch_sample_ids, by:0)
			.ifEmpty {
				log.error "fastp output (params.OVERRIDE_LIST_FASTP) does not match sample IDs in params.metadata, or row in inputs not spaced with tabs"
				fastp_ch.view{log.error "Debug fastp_ch: ${it}"}
				ch_sample_ids.view{log.error "Debug ch_sample_ids: ${it}"}
				System.exit(1)
			}

			preprocessed_reads = preprocessed_reads.combine(ch_sample_ids, by:0)
				.ifEmpty {
					log.error "preprocessed_reads (params.OVERRIDE_PREPROCESSED) does not match sample IDs in params.metadata, or row in inputs not spaced with tabs"
					preprocessed_reads.view{log.error "Debug preprocessed_reads: ${it}"}
					ch_sample_ids.view{log.error "Debug ch_sample_ids: ${it}"}
					System.exit(1)
				}

			kraken_out = kraken_out.combine(ch_sample_ids, by:0)
				.ifEmpty {
					log.error "kraken_out (params.OVERRIDE_LIST_KRAKEN) does not match sample IDs in params.metadata, or row in inputs not spaced with tabs"
					kraken_out.view{log.error "Debug kraken_out: ${it}"}
					ch_sample_ids.view{log.error "Debug ch_sample_ids: ${it}"}
					System.exit(1)
				}

			mapped_bam = mapped_bam.combine(ch_sample_ids, by:0)
				.ifEmpty {
					log.error "mapped_bam (params.OVERRIDE_LIST_BAM) does not match sample IDs in params.metadata, or row in inputs not spaced with tabs"
					mapped_bam.view{log.error "Debug mapped_bam: ${it}"}
					ch_sample_ids.view{log.error "Debug ch_sample_ids: ${it}"}
					System.exit(1)
				}

			bamdam_bam_lca = bamdam_bam_lca.combine( ch_sample_ids ,by:0 )
			.ifEmpty {
				log.error "bamdam_bam_lca (params.OVERRIDE_LIST_BAMDAM) does not match with the sample IDs in params.metadata, or row in inputs not spaced with tabs"
				bamdam_bam_lca.view{log.error "Debug bamdam_bam_lca: ${it}"}
				ch_sample_ids.view{log.error "Debug ch_sample_ids: ${it}"}
				System.exit(1)
			}

			paired_reads
			.combine( fastp_ch ,by:0 )
			.combine( preprocessed_reads ,by:0 )
			.combine( kraken_out_metrics, by:0 )
			.combine( mapped_bam, by:0 )
			.combine( bamdam_bam_lca.map{it -> tuple(it[0], it[1])}, by:0 )
			.set{ input_metrics }
			//input_metrics.view()

			METRICS( input_metrics )
			//METRICS.out.collect().view()
			CONCAT_METRICS( METRICS.out.collect() )

			CONCAT_METRICS.out.set{ metrics1 }

		} else {metrics1 = Channel.empty()}

		if (params.OVERRIDE_LIST_METRICS) {
			Channel.fromPath(params.OVERRIDE_LIST_METRICS)
			.set{ metrics2 }
		} else {metrics2 = Channel.empty()}

		metrics = metrics1.concat(metrics2).unique().ifEmpty {
			ch_sample_ids.map { id -> [id, params.METAJAM_DIR+"/assets/NO_FILE8"] }
		}

		if (params.ENABLE_PLOTS == "enable") {

			// metrics.view{log.info "Debug metrics files for plots: ${it}"}
			// bamdam_bam_lca.map{it[3]}.collect().view{log.info "Debug bamdam_bam_lca files for plots: ${it}"}
			// mmseq2_evaluation.map{it[1]}.collect().view{log.info "Debug mmseq2_evaluation files for plots: ${it}"}

			metrics
			.combine(bamdam_bam_lca.map{it[3]}.collect())
			.combine(mmseq2_evaluation.map{it[1]}.collect())
			.map{it -> tuple(it[0], it[1], it[2],
				params.metadata,
				params.PLOTS_SAMPLES_FOR_PLOTS,
				params.PLOTS_MIN_READS,
				params.PLOTS_MODE,
				params.PLOTS_DAMAGE_THRESHOLD,
				params.PLOTS_PLOT_LOW_DAMAGE_TAXA,
				params.PLOTS_EXCLUDE_TAXA,
				params.PLOTS_TAXA_PER_PLOT,
				params.PLOTS_LIST_TAXA_EVOLUTION_FILE,
				params.MAP_LAST_DB_TAG
			)}
			.set{ input_plots }

			//input_plots.view{"Debug input_plots: ${it}"}

			PLOTS( input_plots )
			//PLOTS_KRONA_BY_SITE seperated from PLOTS for a cleaner conda env
			PLOTS_KRONA_BY_SITE(
				bamdam_bam_lca.map{it -> it[3]}.collect(), 
				params.metadata, 
				params.BAMDAM_MINREADS,
				params.BAMDAM_MAXDAMAGE
			)
		}
	}
}
