#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Define absolute paths to pseudo-reads and annotation
workflow MASK_MICROBIAL_LIKE_REGION {
  take:
    input_dir 
    input_list
    mapping_indexes
    pseudo_reads_file_dir
    type_of_pseudo_reads 
    n_allowed_multimappers
    work_dir
    fna2name
    enable_mask_fasta
  
  main:

    // pseudo_reads_file = Channel.fromPath("${pseudo_reads_file_dir}/*", glob: true)
    pseudo_reads_file = Channel
    .fromPath(pseudo_reads_file_dir, type: 'dir')
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
    if (!mapping_indexes) {

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
      
      input_for_align = index_reference.out
            .combine(pseudo_reads_file)
            .map{it -> tuple(it, type_of_pseudo_reads, n_allowed_multimappers)}
                  .map{it -> it.flatten()}
                  .map{it -> tuple(it[0], tuple(it[1],it[2],it[3],it[4],it[5],it[6]),it[7], type_of_pseudo_reads, n_allowed_multimappers)}

      //input_for_align.view()
    } else {
      input_for_align = mapping_indexes
      .combine(pseudo_reads_file)
      .map{it -> tuple(it, type_of_pseudo_reads, n_allowed_multimappers)}
            .map{it -> it.flatten()}
            .map{it -> tuple(it[0], tuple(it[1],it[2],it[3],it[4],it[5],it[6]),it[7], type_of_pseudo_reads, n_allowed_multimappers)}
    }

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

    bowtie2 --large-index -f -k ${n_allowed_multimappers} -x \${ref_name} \
        --end-to-end --quiet --threads "${task.cpus}" --very-sensitive \
        -U ${input_pseudo_reads} | \
        samtools view -bS -F 4 -h -@ "${task.cpus}" - | \
        samtools sort -@ "${task.cpus}" - > PseudoReads_aligned_to_\${input_pseudo_reads_name}.bam
    
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

        samtools sort -@ "${task.cpus}" -o ${ID}_merged.sorted.bam ${ID}.merged.bam

        samtools index -c ${ID}_merged.sorted.bam
    """
}


// Process 3: Detection
process detect_exogenous {

  publishDir "${params.OUTPUT_Dir}/06_masked_bam", 
        mode: "copy"

  container 'docker://quay.io/biocontainers/mulled-v2-0697a5880de9863c66cba89c8310687052a940fc:c72ea422cf70582757ae5648f79b19857320259b-0'

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
  cut -f 2,3,4 "$raw_bed" | tail -n +2 > "\$out"
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
