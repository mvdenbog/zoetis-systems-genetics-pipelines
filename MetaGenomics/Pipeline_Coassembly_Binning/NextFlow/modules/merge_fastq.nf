process MERGE_FASTQ {
    tag "$meta.id"
    publishDir "${params.outdir}/02_assembly/merged_fastq", mode: 'copy'
    label 'MERGE_FASTQ' 

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    script:
    def readList = reads.collect{ it.toString() }
    def read1 = []
    def read2 = []
    readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
    """
    zcat ${read1.join(' ')} > ${meta.sample}_1.merged.fastq
    zcat ${read2.join(' ')} > ${meta.sample}_2.merged.fastq
    gzip ${meta.sample}_1.merged.fastq 
    gzip ${meta.sample}_2.merged.fastq
    """
}