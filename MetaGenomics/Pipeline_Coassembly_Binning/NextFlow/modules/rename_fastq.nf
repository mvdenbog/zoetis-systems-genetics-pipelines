process RENAME_FASTQ {
    tag "$meta.id"
    label 'RENAME_FASTQ' 

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fastq.gz"), emit: reads

    script:
    """
    if [[ "${meta.type}" == "SR" && "${reads[0]}" != "${meta.id}_R1.fastq.gz" ]]
    then
        mv ${reads[0]} ${meta.id}_R1.fastq.gz
        mv ${reads[1]} ${meta.id}_R2.fastq.gz
    fi
    if [[ "${meta.type}" == "HIFI" && "${reads}" != "${meta.id}.fastq.gz" ]]
    then
        mv ${reads} ${meta.id}.fastq.gz
    fi
    """
}