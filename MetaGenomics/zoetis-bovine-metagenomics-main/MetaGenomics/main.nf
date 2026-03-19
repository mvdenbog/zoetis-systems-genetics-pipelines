#!/usr/bin/env nextflow

/*--------------------------------------------------
  Author: Mathias Vandenbogaert, December 2023
---------------------------------------------------*/


nextflow.enable.dsl = 2

include { DATABASES as PHASE_01_DATABASES  } from './subworkflows/PHASE_01_databases'
include { CLEAN_QC as PHASE_01_CLEAN_QC } from './subworkflows/PHASE_01_clean_QC'
include { GET_SOFTWARE_VERSIONS } from './modules/get_software_versions'
include { MULTIQC } from './modules/multiqc'
include { ASSEMBLY as PHASE_02_ASSEMBLY } from './subworkflows/PHASE_02_ASSEMBLY'

workflow {
  assembly_tool = params.assembly

    if ( assembly_tool == null || assembly_tool == ''){
      assembly_tool = 'metaspades'

    }
    else if (!['metaspades','megahit'].contains(assembly_tool)){
      exit 1, "Invalid short read assembler option: ${assembly_tool}. Valid options for short reads: 'metaspades', 'megahit'"
    }

    if (!(params.skip_host_filter) && !(params.host_fasta) && !(params.skip_clean)) {
      exit 1, "You must specify --host_fasta or skip cleaning host step with option --skip_host_filter or skip all clean and qc modules with --skip_clean"
    }




Channel.fromFilePairs(params.input)
.map { sample,  reads ->
        meta = [
            sample:sample,
        ]
        [ "sample":sample,
          "flowcell":null,
          "group":null,
          "fastq_1":reads.first(),
          "fastq_2":reads.last(),
          "assembly":null
        ]
    }
.set { ch_inputs}


  ch_inputs
    .map { item -> 
        def meta = [:] 
        meta.id = item.sample
        if (item.flowcell!=null) { meta.id = meta.id+"_"+item.flowcell}
        meta.sample = item.sample
        meta.flowcell = item.flowcell
        meta.group = item.group
        meta.assembly = item.assembly!=null
          return [meta,[item.fastq_1,item.fastq_2]]
      }
    .set{ch_reads}


  ch_inputs
    .map { item ->
      def meta = [:] 
      meta.id = item.sample
      if (item.flowcell!=null) { meta.id = meta.id+"_"+item.flowcell}
      meta.sample = item.sample
      meta.flowcell = item.flowcell
      meta.group = item.group
      meta.assembly = item.assembly!=null
      return [meta,item.assembly]
    }
    .set { ch_assembly } 
  has_assembly = null // (file(params.input).splitCsv ( header:true, sep:',' ).assembly[0] != null)
  has_flowcell = null //(file(params.input).splitCsv ( header:true, sep:',' ).flowcell[0] != null)
  
  // Databases
  ch_host_fasta = Channel.empty()
  ch_host_index = Channel.empty()
  ch_kaiju_db = Channel.empty()
  ch_kraken2_db = Channel.empty()
  ch_k2_taxonomy_db = Channel.empty()
  ch_taxonomy = Channel.empty()

  PHASE_01_DATABASES ()
  ch_host_fasta = PHASE_01_DATABASES.out.host_fasta
  ch_host_index = PHASE_01_DATABASES.out.host_index
  ch_kaiju_db = PHASE_01_DATABASES.out.kaiju_db
  ch_kraken2_db = PHASE_01_DATABASES.out.kraken2_db
  ch_k2_taxonomy_db = PHASE_01_DATABASES.out.k2_taxonomy_db
  ch_taxonomy = PHASE_01_DATABASES.out.taxonomy

  ch_multiqc_config = Channel.empty()

  /////////
  ///report
  /////////

  ch_cutadapt_report = Channel.empty()
  ch_sickle_report = Channel.empty()
  ch_before_host_filter_report = Channel.empty()
  ch_fastqc_raw_report = Channel.empty()
  ch_fastqc_clean_report = Channel.empty()
  ch_kaiju_report = Channel.empty()
  ch_kraken2_report = Channel.empty()
  ch_unfilter_assembly_flagstat_report = Channel.empty()
  ch_software_versions = Channel.empty()
  ch_software_versions_PH01 = Channel.empty()
  ch_software_versions_PH02 = Channel.empty()
  ch_software_versions_total = Channel.empty()

  ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)

  if ( !params.skip_clean ) {
    PHASE_01_CLEAN_QC (
      ch_reads,
      ch_host_fasta,
      ch_host_index,
      ch_kaiju_db,
      ch_kraken2_db,
      ch_k2_taxonomy_db
    )
    ch_reads = PHASE_01_CLEAN_QC.out.preprocessed_reads
    ch_cutadapt_report = PHASE_01_CLEAN_QC.out.cutadapt_report
    ch_sickle_report = PHASE_01_CLEAN_QC.out.sickle_report

    ch_before_host_filter_report = PHASE_01_CLEAN_QC.out.before_filter_report

    ch_fastqc_raw_report = PHASE_01_CLEAN_QC.out.fastqc_raw_report
    ch_fastqc_clean_report = PHASE_01_CLEAN_QC.out.fastqc_clean_report
    ch_kaiju_report = PHASE_01_CLEAN_QC.out.kaiju_report
    ch_kaiju_v = PHASE_01_CLEAN_QC.out.v_kaiju
    ch_kronatool_v = PHASE_01_CLEAN_QC.out.v_kronatools
    ch_kaiju_fq_species = PHASE_01_CLEAN_QC.out.kaiju_fq_species

    ch_kraken2_report = PHASE_01_CLEAN_QC.out.kraken2_report

    ch_software_versions_PH01 = PHASE_01_CLEAN_QC.out.software_versions
  }



PHASE_02_ASSEMBLY(ch_kaiju_fq_species, ch_kaiju_db)

ch_software_versions_PH02 = PHASE_02_ASSEMBLY.out.software_versions

  ch_software_versions_total = ch_software_versions_PH01.mix(ch_software_versions_PH02).unique { it.getBaseName() }.collect()


  GET_SOFTWARE_VERSIONS(ch_software_versions_total) 
  ch_software_versions = GET_SOFTWARE_VERSIONS.out.yaml

 MULTIQC (
    ch_multiqc_config,
    ch_software_versions,
    ch_cutadapt_report.collect().ifEmpty([]),
    ch_sickle_report.collect().ifEmpty([]),
    ch_before_host_filter_report.collect().ifEmpty([]),
    ch_fastqc_raw_report.collect().ifEmpty([]),
    ch_fastqc_clean_report.collect().ifEmpty([]),
    ch_kaiju_report.collect().ifEmpty([]),
    ch_kraken2_report.collect{it[1]}.ifEmpty([]),
    ch_unfilter_assembly_flagstat_report.collect().ifEmpty([])

  )
  multiqc_report = MULTIQC.out.report

}
