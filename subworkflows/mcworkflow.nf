#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define absolute paths to pseudo-reads and annotation
workflow MASK_MICROBIAL_LIKE_REGION {
  take:
    input_dir 
    input_list
    BOWTIE2_MAPPING_DB1
    BOWTIE2_MAPPING_DB2
    BOWTIE2_MAPPING_DB3
    BOWTIE2_MAPPING_DB4
    BOWTIE2_MAPPING_DB5
    BOWTIE2_MAPPING_DB6
    BOWTIE2_MAPPING_DB7
    BOWTIE2_MAPPING_DB8
    BOWTIE2_MAPPING_DB9
    BOWTIE2_MAPPING_DB10
    pseudo_reads_file_dir
    type_of_pseudo_reads 
    n_allowed_multimappers
    work_dir
    fna2name
    enable_mask_fasta
  
  main:

    pseudo_reads_file = Channel.fromPath(pseudo_reads_file_dir, type: 'dir')
        .flatMap { dir -> dir.listFiles() }
        .filter {
            it.name.endsWith('.fna') ||
            it.name.endsWith('.fa') ||
            it.name.endsWith('.fasta') ||
            it.name.endsWith('.fna.gz') ||
            it.name.endsWith('.fa.gz') ||
            it.name.endsWith('.fasta.gz')
        }
    .ifEmpty {
        log.error "No input files found in ${pseudo_reads_file_dir}"
        System.exit(1)
    }

    // prioritize reading mapping indexes if provided, otherwise generate them from input fasta files
    if (!params.BOWTIE2_MAPPING_DB1) {

      input1 = Channel.empty()
      if (input_dir != ""){
        input1 = Channel
        .fromPath(input_dir, type: 'dir')
        .flatMap { dir -> dir.listFiles() }
        .filter { 
            it.name.endsWith('.fna') ||
            it.name.endsWith('.fa') ||
            it.name.endsWith('.fasta') ||
            it.name.endsWith('.fna.gz') ||
            it.name.endsWith('.fa.gz') ||
            it.name.endsWith('.fasta.gz')
        }
        .ifEmpty {
          log.error "No input files (.fna, .fa, .fasta, optionally .gz) found in ${input_dir}"
          System.exit(1)
        }
        .map { f ->tuple(f.baseName.replaceFirst(/(\.fna|\.fa|\.fasta)(\.gz)?$/, ''),f)}
        // .view()
      }
      
      input2 = Channel.empty()
      if (input_list != ""){
        input2 = Channel
        .fromPath(input_list)
        .splitText()
        .map { it.trim() }
        .filter { it }
        .map { f ->file(f)}
        .ifEmpty {
          log.error "No input files (.fna, .fa, .fasta, optionally .gz) found in ${input_list}"
          System.exit(1)
        }
        .map { f ->tuple(f.baseName.replaceFirst(/(\.fna|\.fa|\.fasta)(\.gz)?$/, ''),f)}
        // .view()
      }
      
      input1.concat(input2).unique().set{files}
      //files.view()

      index_reference(files)

      index_reference.out.set{mapping_indexes}

    } else {
      // Channel.fromPath( BOWTIE2_MAPPING_DBs )
			// .splitText { it.strip( ) }
			// .map { it -> 
			// def name = it.tokenize('/')[-1]   // get basename
			// tuple(name, file("${it}*bt2*"))}
			// .groupTuple()
			// .map { idx, idxs -> tuple(idx, idxs[0]) }
			// .set { mapping_indexes }

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
    }
    
    input_for_align = mapping_indexes
      .combine(pseudo_reads_file)
      .map{it -> tuple(it, type_of_pseudo_reads, n_allowed_multimappers)}
            .map{it -> it.flatten()}
            .map{it -> tuple(it[0], tuple(it[1],it[2],it[3],it[4],it[5],it[6]),it[7], type_of_pseudo_reads, n_allowed_multimappers)}


    //input_for_align.view{"Input for alignment: ${it}"}

    align_pseudo_reads(input_for_align)

    merge_bam(align_pseudo_reads.out.groupTuple())

    merge_bam.out
    .map{it -> tuple(it, type_of_pseudo_reads, work_dir, fna2name)}
    .map{it -> it.flatten()}
    .set{  detect_input }

    // detect_exogenous(align_pseudo_reads.out.bam, align_pseudo_reads.out.ref, type_of_pseudo_reads, work_dir, fna2name)
    detect_exogenous(detect_input)

    make_bedfile(detect_exogenous.out.for_bedfile)

    
    if (enable_mask_fasta == "enable" ) {

      if (!input_dir && !input_list) {
        log.error "Error: To enable masking fasta, you must provide input fasta files through either params.MCWORKFLOW_input_dir or params.MCWORKFLOW_input_list."
        System.exit(1)
      } else {
        //make_bedfile.out.combine( files, by:0 ).view()
        mask_fasta( make_bedfile.out.combine( files,by:0 ) )
      }
    } 

  emit:
    make_bedfile.out

}

// Process 1: Indexing
process index_reference {
    conda './envs/bowtie2.yml'

    input:
    tuple val(ID), path(input_ref)

    output:
    tuple val(ID), path("*.bt2l")

    script:
    """
    bowtie2-build --large-index \$(basename ${input_ref}) \$(basename ${input_ref}) --threads "${task.cpus}"
    """
}

// Process 2: Alignment
process align_pseudo_reads {
    conda './envs/bowtie2.yml'

    input:
    tuple val(ID), path(index), path(input_pseudo_reads), val(type_of_pseudo_reads), val(n_allowed_multimappers)

    output:
    tuple val(ID), path("*.bam")

    script:
    """
    index1=\$(printf '%s\n' *.bt2* | head -n1)
    ref_name=\$(echo \$index1 | sed 's/.1.bt2l//' | sed 's/.1.bt2//')

    input_pseudo_reads_name=\$(basename "$input_pseudo_reads" | sed -E 's/\\.(fna|fa|fasta)(\\.gz)?\$//')

    if bowtie2 --large-index -f -k ${n_allowed_multimappers} -x \${ref_name} \
        --end-to-end --quiet --threads "${task.cpus}" --very-sensitive \
        -U ${input_pseudo_reads} | \
        samtools view -bS -F 4 -h -@ "${task.cpus}" - | \
        samtools sort -@ "${task.cpus}" - > PseudoReads_aligned_to_\${input_pseudo_reads_name}.bam; then
		echo "map successfully with large index"
    else
	#in case of small index .bt2
	bowtie2 -f -k ${n_allowed_multimappers} -x \${ref_name} \
        --end-to-end --quiet --threads "${task.cpus}" --very-sensitive \
        -U ${input_pseudo_reads} | \
        samtools view -bS -F 4 -h -@ "${task.cpus}" - | \
        samtools sort -@ "${task.cpus}" - > PseudoReads_aligned_to_\${input_pseudo_reads_name}.bam
	echo "map successfully with small index"
    fi
    
    input_ref_name=\$(echo \$index1 | sed 's/.1.bt2l//' | sed -E 's/.fasta|.fa|.fna//' | sed 's/.gz//')
    """
}


// merge bam of same sample mapped to different databases
process merge_bam {
    
    conda './envs/bowtie2.yml'
    input:
        tuple val(ID), path(bams)

    output:
        tuple val(ID), path("*_merged.sorted.bam"), path("*.bam.csi")

    script:
    """
        mapped_bams=()        
        for bam1 in *.bam; do
            [[ "\$bam1" == *mapped.bam ]] && continue
            out_bam="\$(basename \$bam1 .bam).mapped.bam"
            samtools view -@ 80 -b -F 0x4 "\$bam1" -o "\$out_bam"
            if [[ \$(samtools view -c "\$out_bam") -ne 0 ]]; then mapped_bams+=("\$out_bam"); fi   # only keep non-empty BAMs
        done

        if [ \${#mapped_bams[@]} -eq 0 ]; then
            echo "ERROR: No mapped reads found in any BAM." >&2
            exit 1
        fi

        samtools merge ${ID}.merged.bam "\${mapped_bams[@]}"

        samtools quickcheck ${ID}.merged.bam || {
            echo "ERROR: Merging is not successful: ${ID}.merged.bam" >&2
            exit 1
        }

        samtools sort -@ "${task.cpus}" -o ${ID}_merged.sorted.bam ${ID}.merged.bam

        samtools index -c ${ID}_merged.sorted.bam
    """
}


// Process 3: Detection
process detect_exogenous {

  publishDir "${params.OUTPUT_Dir}/06_masked_bam", 
        mode: "copy"

  conda './envs/samtools_R.yml'

  input:
    tuple val(input_ref), path(bam), path(bai), val(type_of_pseudo_reads), val(work_dir), path(fna2name)

  output:
    path("*abund_*.txt")
    tuple val(input_ref), path("*coords_micr_like_regions*.txt"), emit: for_bedfile
    path("*boc_*.txt")
    path("*_microbes_abundant_*.txt")
    
    

  script:
  """
  #get just bam
  bamfile=\$(echo $bam | awk '{print \$1}' )

  detect_exogenous.sh \
      \${bamfile} \
      ${input_ref} \
      ${type_of_pseudo_reads} \
      ${work_dir}

  echo "GENERATE COORDINATIONS OF MICROBIAL-LIKE REGIONS (BEDFILES)"
  for j in \$(cat refs_uniq_sorted.txt)
	do
	echo \${j} CONTIG OF ${input_ref}
	extract_coords.R ${type_of_pseudo_reads} \${j}__${input_ref}.boc $fna2name
	echo DELETING BAM AND COMPRESSING BOC FILES
	rm \${j}.bam
	rm \${j}__${input_ref}.boc
  done
  #remove intermediate files
  rm refs_uniq_sorted.txt refs_uniq_sorted_reads.txt total_length_per_ref.txt boc_per_ref.txt
  rm $fna2name #avoid output it

  #add prefix
  for f in \$(ls * | grep -v .bam); do
	[[ "\$f" == *.bam ]] && continue
    mv "\$f" "${input_ref}_${type_of_pseudo_reads}_\$f"
  done

  """
}

process make_bedfile {

  publishDir "${params.OUTPUT_Dir}/06_masked_bam", mode: "copy"

  input: 
    tuple val(ID), path(raw_bed)

  output: 
    tuple val(ID), path("*.bed")

  script:
  """
  out=\$(basename "$raw_bed" .txt).bed
  cut -f 2,3,4 "$raw_bed" | tail -n +2 | awk '{
    \$2=sprintf("%.0f",\$2);
    \$3=sprintf("%.0f",\$3);
    print}' OFS='\t' > "\$out"
  """
}

// mask the fasta with bedfile
process mask_fasta {

  conda 'bioconda::bedtools'

  publishDir "${params.OUTPUT_Dir}/06_masked_bam", mode: "copy"
  
  input: 
  tuple val(ID), path(bed), path(ref)

  output: path("*.masked.fna")

  script:
  """
    bedtools maskfasta -fi ${ref} -bed ${bed} -fo ${ID}.masked.fna
  """
}
