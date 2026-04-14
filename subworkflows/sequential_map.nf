workflow {

    take:
        mapping_indexes_raw
        reads

    main:
        mapping_indexes_raw
            .toList()
            .set { mapping_indexes }

        reads
            .map { it -> tuple(it[0], it[1], "10") }
            .combine(mapping_indexes.map{it -> it[0]})
            .set{input_bowtie2}

        // collector
        all_outputs = Channel.empty()

        BOWTIE2_DB1(input_bowtie2)
        all_outputs = all_outputs.mix(BOWTIE2_DB1.out.mapped.map{it -> it[1]})

        if (mapping_indexes.map{it -> it[1]}) {

            BOWTIE2_DB1.out.unmapped
                .map { it -> tuple(it[0], it[1], "10") }
                .combine(mapping_indexes.map{it -> it[1]})
                .set{input_bowtie2_db2}

            BOWTIE2_DB2(input_bowtie2_db2)
            all_outputs = all_outputs.mix(BOWTIE2_DB2.out.mapped.map{it -> it[1]})

            if (mapping_indexes.map{it -> it[2]}) {

                BOWTIE2_DB2.out.unmapped
                    .map { it -> tuple(it[0], it[1], "10") }
                    .combine(mapping_indexes.map{it -> it[2]})
                    .set{input_bowtie2_db3}

                BOWTIE2_DB3(input_bowtie2_db3)
                all_outputs = all_outputs.mix(BOWTIE2_DB3.out.mapped.map{it -> it[1]})

                if (mapping_indexes.map{it -> it[3]}) {

                    BOWTIE2_DB3.out.unmapped
                        .map { it -> tuple(it[0], it[1], "10") }
                        .combine(mapping_indexes.map{it -> it[3]})
                        .set{input_bowtie2_db4}

                    BOWTIE2_DB4(input_bowtie2_db4)
                    all_outputs = all_outputs.mix(BOWTIE2_DB4.out.mapped.map{it -> it[1]})

                    if (mapping_indexes.map{it -> it[4]}) {

                        BOWTIE2_DB4.out.unmapped
                            .map { it -> tuple(it[0], it[1], "10") }
                            .combine(mapping_indexes.map{it -> it[4]})
                            .set{input_bowtie2_db5}

                        BOWTIE2_DB5(input_bowtie2_db5)
                        all_outputs = all_outputs.mix(BOWTIE2_DB5.out.mapped.map{it -> it[1]})

                        if (mapping_indexes.map{it -> it[5]}) {

                            BOWTIE2_DB5.out.unmapped
                                .map { it -> tuple(it[0], it[1], "10") }
                                .combine(mapping_indexes.map{it -> it[5]})
                                .set{input_bowtie2_db6}

                            BOWTIE2_DB6(input_bowtie2_db6)
                            all_outputs = all_outputs.mix(BOWTIE2_DB6.out.mapped.map{it -> it[1]})

                            if (mapping_indexes.map{it -> it[6]}) {

                                BOWTIE2_DB6.out.unmapped
                                    .map { it -> tuple(it[0], it[1], "10") }
                                    .combine(mapping_indexes.map{it -> it[6]})
                                    .set{input_bowtie2_db7}

                                BOWTIE2_DB7(input_bowtie2_db7)
                                all_outputs = all_outputs.mix(BOWTIE2_DB7.out.mapped.map{it -> it[1]})

                                if (mapping_indexes.map{it -> it[7]}) {

                                    BOWTIE2_DB7.out.unmapped
                                        .map { it -> tuple(it[0], it[1], "10") }
                                        .combine(mapping_indexes.map{it -> it[7]})
                                        .set{input_bowtie2_db8}

                                    BOWTIE2_DB8(input_bowtie2_db8)
                                    all_outputs = all_outputs.mix(BOWTIE2_DB8.out.mapped.map{it -> it[1]})

                                    if (mapping_indexes.map{it -> it[8]}) {

                                        BOWTIE2_DB8.out.unmapped
                                            .map { it -> tuple(it[0], it[1], "10") }
                                            .combine(mapping_indexes.map{it -> it[8]})
                                            .set{input_bowtie2_db9}

                                        BOWTIE2_DB9(input_bowtie2_db9)
                                        all_outputs = all_outputs.mix(BOWTIE2_DB9.out.mapped.map{it -> it[1]})

                                        if (mapping_indexes.map{it -> it[9]}) {

                                            BOWTIE2_DB9.out.unmapped
                                                .map { it -> tuple(it[0], it[1], "10") }
                                                .combine(mapping_indexes.map{it -> it[9]})
                                                .set{input_bowtie2_db10}

                                            BOWTIE2_DB10(input_bowtie2_db10)
                                            all_outputs = all_outputs.mix(BOWTIE2_DB10.out.mapped.map{it -> it[1]})
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }

    all_outputs.view()

    output:
        all_outputs.flatten().collect()
}


process BOWTIE2_DB1 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)
    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper
    """
}


process BOWTIE2_DB2 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB3 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB4 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB5 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB6 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB7 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB8 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB9 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB10 {
    // label 'small_memory' // for test
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}
