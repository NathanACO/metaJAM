#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// include the subworkflow
include { FASTP; SGA; PRINSEQ; KRAKEN2; BOWTIE2; MERGE_BAM;  MASK_REGIONS; FILTERBAM; NGSLCA;  BAMDAM; MMSEQ2 ; KRONATOOLS; METRICS; CONCAT_METRICS; PLOTS } from './func.nf'

workflow {
	// input =================
	Channel.fromFilePairs(params.FASTQ)
	.map { sample_id, reads -> tuple(sample_id, reads[0], reads[1])}
	.set { paired_reads }

	Channel.fromPath( params.BOWTIE2_MAPPING_DBs )
	.splitText { it.strip( ) }
	.map { it -> 
    def name = it.tokenize('/')[-1]   // get basename
    tuple(name, file("${it}*bt2*"))}
	.groupTuple()
	.map { idx, idxs -> tuple(idx, idxs[0]) }
	.set { mapping_indexes }
	// mapping_indexes.view()

	def acc2taxid = file(params.ACC2TAXID) 
	// def acc2taxid_exists = acc2taxid.exists()

	//add check if file exist or empty

	input = paired_reads.map { tuple(it[0], it[1], it[2], params.FASTP_OVERLAP_LEN_REQUIRE, params.FASTP_MIN_LENGTH) }

	if (params.ENABLE_PREPROCESS == "enable") {
		if (params.ENABLE_FASTP == "enable") { FASTP( input ) }

    	if (params.ENABLE_PRINSEQ == "enable") {  
			SGA(FASTP.out .map { tuple(it[0], it[1], params.SGA_DUST_THRESHOLD) })
			// SGA.out.rm_low_complexity.view()
			 }
		if (params.ENABLE_PRINSEQ == "enable") {  
			PRINSEQ(FASTP.out .map { tuple(it[0], it[1], params.PRINSEQ_COMPLEXITY_METHOD, params.PRINSEQ_COMPLEXITY_THRESHOLD, params.PRINSEQ_MIN_LEN, params.PRINSEQ_DEREP)}) 
			// PRINSEQ.out.rm_low_complexity.view()
			}
		if (params.ENABLE_LOW_COMPLEXITY_FILTER=="SGA"){
			SGA.out.rm_low_complexity.set { preprocessed_reads }
		} else if (params.ENABLE_LOW_COMPLEXITY_FILTER=="PRINSEQ") {
			PRINSEQ.out.rm_low_complexity.set { preprocessed_reads }
		}
	}

	// preprocessed_reads.view()

	if (params.ENABLE_KRAKEN_GTDB == "enable") { KRAKEN2( preprocessed_reads.map { tuple(it[0], it[1], params.KRAKEN2_FILTER_DATABASE) }) }
	// KRAKEN2.out.not_microbe.view()

	KRAKEN2.out.not_microbe
	.map { it -> tuple(it[0], it[1], params.BOWTIE2_N_ALLOW_MULTIMAPPER) }
	.combine( mapping_indexes ) 
	.set { input_mapping }
	
	// input_mapping.view()

	if (params.ENABLE_MAPPING == "enable") { BOWTIE2( input_mapping )}
	// BOWTIE2.out.groupTuple().view() //group by sample ID

	acc2taxid = file( params.ACC2TAXID ) 

	//BOWTIE2.out.groupTuple().view()

	//prioritize processing mapped_bam if specified than bowtie2 mapping results
	if ( params.MAPPED_BAM == "" ) {
		BOWTIE2.out.groupTuple().set { mapped_bam }
	} else {
		Channel.fromPath( params.MAPPED_BAM , checkIfExists: true)
		.splitCsv( sep: "\t" )
		.map { row -> [ row[0], file(row[1]) ] }
		.groupTuple()
		.set { mapped_bam }
	}
	//mapped_bam.view()

	//generate acc2taxid if not provided // to test
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

	// MERGE_BAM( BOWTIE2.out.groupTuple() )
	MERGE_BAM( mapped_bam )
	// MERGE_BAM.out.view()

	// MERGE_BAM.out.map { tuple(it[0], it[1], params.REGIONS_TO_MASK) } .view()

	if (params.ENABLE_MASK_REGIONS == "enable") { 
		MASK_REGIONS( MERGE_BAM.out.map { tuple(it[0], it[1], params.REGIONS_TO_MASK) }  )
		MASK_REGIONS.out.set{ bam }
	} else {
		MERGE_BAM.out.set{ bam }
	}
	
	// if (params.ENABLE_FILTERBAM == "enable") { 
	// 	FILTERBAM( input )
	// 	FILTERBAM.out.set{ bam }
	// 	} else { MERGE_BAM.out.set{ bam } }

	// bam.map { tuple(it[0], it[1], params.NAMES, params.NODES, params.ACC2TAXID )}.view()

	if (params.ENABLE_NGSLCA == "enable") { NGSLCA( bam.map { tuple(it[0], it[1], params.NAMES, params.NODES, acc2taxid )} )}
	// NGSLCA.out.lca.view()

	MERGE_BAM.out
	.combine( NGSLCA.out.lca ,by:0 )
	.map(it -> tuple(it[0], it[1], it[2],
	params.BAMDAM_STRANDED,
	params.BAMDAM_MINREADS,
	params.BAMDAM_MAXDAMAGE,
	params.BAMDAM_TOP_GENUS))
	.set{ input_bamdam }

	//input_bamdam.view()
	BAMDAM( input_bamdam )

	KRONATOOLS( BAMDAM.out.xml )

	MMSEQS2_GENERA_FILE = Channel.fromPath(params.MMSEQS2_GENERA_FILE, checkIfExists:true)

	BAMDAM.out.lca
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

	BAMDAM.out.lca.map{it -> tuple(it[0], it[1])}.set{bamdam_bam}

	paired_reads
	.combine( FASTP.out ,by:0 )
	.combine( preprocessed_reads ,by:0 )
	.combine( KRAKEN2.out.not_microbe ,by:0 )
	.combine( BOWTIE2.out.groupTuple(), by:0 )
	.combine( bamdam_bam, by:0 )
	.set{ input_metrics }
	// input_metrics.view()

	METRICS( input_metrics )
	METRICS.out.collect().view()
	CONCAT_METRICS( METRICS.out.collect() )

	CONCAT_METRICS.out
	.combine(BAMDAM.out.lca.map(it -> it[2]).collect())
	.combine(MMSEQ2.out.evaluation.collect())
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

