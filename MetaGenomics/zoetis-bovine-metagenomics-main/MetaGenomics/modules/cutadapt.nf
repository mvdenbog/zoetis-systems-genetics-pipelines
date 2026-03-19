process CUTADAPT {
        label 'GENERAL'
	tag "${meta.id}"

	publishDir "${params.outdir}/${meta.id}/CleanedReads", mode: 'copy', pattern: 'cleaned_*.fastq.gz'
	publishDir "${params.outdir}/${meta.id}/CleanedReads/logs", mode: 'copy', pattern: '*_cutadapt.log'

	input:
	tuple val(meta), path(reads)
	val adapter1
	val adapter2

	output:
	tuple val(meta), path("*${meta.id}*.fastq.gz"), emit: reads
	path "${meta.id}_cutadapt.log", emit: report
	path "v_cutadapt.txt", emit: v_cutadapt 

	script:
	if (!params.use_sickle & params.skip_host_filter) {
		// output are final cleaned paths
		output_paths = "-o cleaned_${meta.id}_R1.fastq.gz -p cleaned_${meta.id}_R2.fastq.gz"
	}
	else {
		// tempory paths not saved in publish dir
		output_paths = "-o ${meta.id}_cutadapt_R1.fastq.gz -p ${meta.id}_cutadapt_R2.fastq.gz"
	}
	if (!params.use_sickle){
		quality_trim =  "-q 20,20" 
	} else {
		quality_trim = ""
	}	
	"""
	cutadapt -a $adapter1 -A $adapter2 $output_paths -m 36 --trim-n -q 20,20 --max-n 0 \
	--cores=${task.cpus} ${reads[0]} ${reads[1]} > ${meta.id}_cutadapt.log

	cutadapt --version &> v_cutadapt.txt
	"""		
}
