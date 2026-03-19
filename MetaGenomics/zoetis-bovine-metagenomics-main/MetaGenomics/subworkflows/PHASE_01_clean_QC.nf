/*--------------------------------------------------
  Author: Mathias Vandenbogaert, December 2023
---------------------------------------------------*/


include { CUTADAPT } from '../modules/cutadapt'
include { SICKLE } from '../modules/sickle'
include { HOST_FILTER } from '../modules/host_filter'
include { FASTQC as FASTQC_RAW; FASTQC as FASTQC_CLEANED } from '../modules/fastqc'
include { KAIJU_AND_MERGE } from '../modules/kaiju'
include { KRAKEN2 } from '../modules/kraken2'
include { KRAKEN2ToKRONA } from '../modules/kraken2'
include { KRAKEN2KRONAMERGE } from '../modules/kraken2'

workflow CLEAN_QC {
	take:
		raw_reads
		host_fasta
		host_index
		kaiju_db
                kraken2_db
                k2_taxonomy_db
                

	main:
		ch_adapter1 = Channel.value(params.adapter1)
		ch_adapter2 = Channel.value(params.adapter2)

		ch_sickle_report = Channel.empty()
		ch_cutadapt_report = Channel.empty()
		ch_before_filter_report = Channel.empty()
		ch_preprocessed_reads = Channel.empty()
		ch_fastqc_clean_report = Channel.empty()
		ch_kaiju_report = Channel.empty()
                ch_kraken2_report = Channel.empty()

		ch_cutadapt_v = Channel.empty()
		ch_sickle_v = Channel.empty()
		ch_bwa_mem_v = Channel.empty()
		ch_samtools_v = Channel.empty()
		ch_sickle_v = Channel.empty()
		ch_minimap_v = Channel.empty()
		ch_kaiju_v = Channel.empty()
                ch_kraken2_v = Channel.empty()
		ch_kronatool_v = Channel.empty()

			CUTADAPT (
				raw_reads,
				ch_adapter1,
				ch_adapter2
			)
			ch_intermediate_reads = CUTADAPT.out.reads
			ch_cutadapt_report = CUTADAPT.out.report

			ch_cutadapt_v =  CUTADAPT.out.v_cutadapt

			if (params.use_sickle) {
				SICKLE (ch_intermediate_reads)
				ch_intermediate_reads = SICKLE.out.reads
				ch_sickle_report = SICKLE.out.report

				ch_sickle_v = SICKLE.out.v_sickle
			}
			ch_preprocessed_reads = ch_intermediate_reads
		

		if (!params.skip_host_filter) {
				HOST_FILTER (
					ch_intermediate_reads,
					host_fasta,
					host_index
				)
				ch_preprocessed_reads = HOST_FILTER.out.reads

				ch_before_filter_report = HOST_FILTER.out.nf_report

				ch_bwa_mem_v = HOST_FILTER.out.v_bwa_mem2
				ch_samtools_v = HOST_FILTER.out.v_samtools
		}

		FASTQC_RAW(raw_reads)
		FASTQC_CLEANED(ch_preprocessed_reads)
		ch_fastqc_raw_report = FASTQC_RAW.out.zip
		ch_fastqc_clean_report = FASTQC_CLEANED.out.zip

		ch_fastqc_v = FASTQC_RAW.out.v_fastqc

		if (params.skip_host_filter && params.type.toUpperCase() == "HIFI") {
			ch_preprocessed_reads = raw_reads
		}

		if (!params.skip_kaiju) {
			KAIJU_AND_MERGE(
				ch_preprocessed_reads,
				kaiju_db
			)
			ch_kaiju_report = KAIJU_AND_MERGE.out.report
			ch_kaiju_v= KAIJU_AND_MERGE.out.v_kaiju
			ch_kronatool_v=KAIJU_AND_MERGE.out.v_kronatools
                        ch_kaiju_fq_species = KAIJU_AND_MERGE.out.kaiju_fq_species
		} 

                if (!params.skip_kraken2){
                      KRAKEN2(ch_preprocessed_reads, kraken2_db)
                      ch_kraken2_report = KRAKEN2.out.report
                      ch_kraken2_v = KRAKEN2.out.v_kraken2
KRAKEN2ToKRONA(ch_kraken2_report, k2_taxonomy_db)
KRAKEN2KRONAMERGE(KRAKEN2ToKRONA.out.krona_in.collect(), k2_taxonomy_db)
                }

		ch_software_versions = ch_cutadapt_v.first().mix(ch_sickle_v.first(),
												ch_bwa_mem_v.first(),
												ch_samtools_v.first(),
												ch_minimap_v.first(),
												ch_fastqc_v.first(),
												ch_kaiju_v.first(),
                                                                                                ch_kraken2_v.first(),
												ch_kronatool_v
												)



	emit:
		preprocessed_reads = ch_preprocessed_reads
		cutadapt_report = ch_cutadapt_report
		sickle_report = ch_sickle_report
		before_filter_report = ch_before_filter_report
		fastqc_raw_report = ch_fastqc_raw_report
		fastqc_clean_report = ch_fastqc_clean_report
		kaiju_report = ch_kaiju_report
                v_kaiju = ch_kaiju_v
                v_kronatools = ch_kronatool_v
                kaiju_fq_species = ch_kaiju_fq_species
                kraken2_report = ch_kraken2_report
		software_versions = ch_software_versions
}
