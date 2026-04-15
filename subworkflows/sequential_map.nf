workflow SEQUENTIAL_MAP {

    take:
        // mapping_indexes_raw
        reads
        n_allow_multimapper

    main:
        

        // mapping_indexes_raw
        //     .toList()
        //     .set { mapping_indexes }

        // mapping databases
        mapping_index1  = Channel.value(params.BOWTIE2_MAPPING_DB1).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) } 
        reads
            .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
            .combine(mapping_index1)
            .set{input_bowtie2}

        // collector
        all_outputs = Channel.empty()
        // placeholders for outputs
        mapped_db1 = Channel.empty()
        mapped_db2 = Channel.empty()
        mapped_db3 = Channel.empty()
        mapped_db4 = Channel.empty()
        mapped_db5 = Channel.empty()
        mapped_db6 = Channel.empty()
        mapped_db7 = Channel.empty()
        mapped_db8 = Channel.empty()
        mapped_db9 = Channel.empty()
        mapped_db10 = Channel.empty()

        BOWTIE2_DB1(input_bowtie2)
        BOWTIE2_DB1.out.mapped.set{mapped_db1}

        if (params.BOWTIE2_MAPPING_DB2) {

            mapping_index2  = Channel.value(params.BOWTIE2_MAPPING_DB2).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
            BOWTIE2_DB1.out.unmapped
                .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                .combine(mapping_index2)
                .set{input_bowtie2_db2}

            BOWTIE2_DB2(input_bowtie2_db2)
            BOWTIE2_DB2.out.mapped.set{mapped_db2}

            if (params.BOWTIE2_MAPPING_DB3) {
                mapping_index3  = Channel.value(params.BOWTIE2_MAPPING_DB3).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
                BOWTIE2_DB2.out.unmapped
                    .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                    .combine(mapping_index3)
                    .set{input_bowtie2_db3}

                BOWTIE2_DB3(input_bowtie2_db3)
                BOWTIE2_DB3.out.mapped.set{mapped_db3}

                if (params.BOWTIE2_MAPPING_DB4) {

                    mapping_index4  = Channel.value(params.BOWTIE2_MAPPING_DB4).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
                    BOWTIE2_DB3.out.unmapped
                        .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                        .combine(mapping_index4)
                        .set{input_bowtie2_db4}

                    BOWTIE2_DB4(input_bowtie2_db4)
                    BOWTIE2_DB4.out.mapped.set{mapped_db4}

                    if (params.BOWTIE2_MAPPING_DB5) {

                        mapping_index5  = Channel.value(params.BOWTIE2_MAPPING_DB5).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
                        BOWTIE2_DB4.out.unmapped
                            .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                            .combine(mapping_index5)
                            .set{input_bowtie2_db5}

                        BOWTIE2_DB5(input_bowtie2_db5)
                        BOWTIE2_DB5.out.mapped.set{mapped_db5}

                        if (params.BOWTIE2_MAPPING_DB6) {

                            mapping_index6  = Channel.value(params.BOWTIE2_MAPPING_DB6).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
                            BOWTIE2_DB5.out.unmapped
                                .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                                .combine(mapping_index6)
                                .set{input_bowtie2_db6}

                            BOWTIE2_DB6(input_bowtie2_db6)
                            BOWTIE2_DB6.out.mapped.set{mapped_db6}

                            if (params.BOWTIE2_MAPPING_DB7) {

                                mapping_index7  = Channel.value(params.BOWTIE2_MAPPING_DB7).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
                                BOWTIE2_DB6.out.unmapped
                                    .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                                    .combine(mapping_index7)
                                    .set{input_bowtie2_db7}

                                BOWTIE2_DB7(input_bowtie2_db7)
                                BOWTIE2_DB7.out.mapped.set{mapped_db7}

                                if (params.BOWTIE2_MAPPING_DB8) {

                                    mapping_index8  = Channel.value(params.BOWTIE2_MAPPING_DB8).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
                                    BOWTIE2_DB7.out.unmapped
                                        .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                                        .combine(mapping_index8)
                                        .set{input_bowtie2_db8}

                                    BOWTIE2_DB8(input_bowtie2_db8)
                                    BOWTIE2_DB8.out.mapped.set{mapped_db8}

                                    if (params.BOWTIE2_MAPPING_DB9) {

                                        mapping_index9  = Channel.value(params.BOWTIE2_MAPPING_DB9).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
                                        BOWTIE2_DB8.out.unmapped
                                            .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                                            .combine(mapping_index9)
                                            .set{input_bowtie2_db9}

                                        BOWTIE2_DB9(input_bowtie2_db9)
                                        BOWTIE2_DB9.out.mapped.set{mapped_db9}

                                        if (params.BOWTIE2_MAPPING_DB10) {

                                            mapping_index10  = Channel.value(params.BOWTIE2_MAPPING_DB10).map { it -> def name = file(it).name; tuple(name, file("${it}*.bt2*"))}.groupTuple().map { idx, idxs -> tuple(idx, idxs[0]) }
                                            BOWTIE2_DB9.out.unmapped
                                                .map { it -> tuple(it[0], it[1], n_allow_multimapper) }
                                                .combine(mapping_index10)
                                                .set{input_bowtie2_db10}

                                            BOWTIE2_DB10(input_bowtie2_db10)
                                            BOWTIE2_DB10.out.mapped.set{mapped_db10}
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

        }
    all_outputs = all_outputs.mix(mapped_db1)
    all_outputs = all_outputs.mix(mapped_db2)
    all_outputs = all_outputs.mix(mapped_db3)
    all_outputs = all_outputs.mix(mapped_db4)
    all_outputs = all_outputs.mix(mapped_db5)
    all_outputs = all_outputs.mix(mapped_db6)
    all_outputs = all_outputs.mix(mapped_db7)
    all_outputs = all_outputs.mix(mapped_db8)
    all_outputs = all_outputs.mix(mapped_db9)
    all_outputs = all_outputs.mix(mapped_db10)

    emit:
        all_outputs
}


process BOWTIE2_DB1 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)
    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper
    """
}


process BOWTIE2_DB2 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB3 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB4 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB5 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB6 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB7 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB8 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB9 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}

process BOWTIE2_DB10 {
    label 'mapping'
    conda './envs/bowtie2.yml'
    input:    
        tuple val(ID), path(reads), val(n_allow_multimapper), 
              val(idx), path(idxs)

    output:
        tuple val(ID), path("*.bam"), emit: mapped
        tuple val(ID), path("${ID}_${idx}_unclass.fq.gz"), emit: unmapped

    publishDir "${params.OUTPUT_Dir}/04_sequantial_mapping", mode: "copy"
        
    script:
    """
    mapping.sh $ID $idx $reads "${task.cpus}" $n_allow_multimapper

    """
}
