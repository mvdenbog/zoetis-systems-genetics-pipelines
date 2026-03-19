#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { DATABASES } from './subworkflows/00_databases'
include { STEP_01_CLEAN_QC as S01_CLEAN_QC } from './subworkflows/01_clean_qc'
include { STEP_02_ASSEMBLY as S02_ASSEMBLY } from './subworkflows/02_assembly'
include { STEP_03_ASSEMBLY_FILTER as S03_FILTERING } from './subworkflows/03_filtering'
include { STEP_04_STRUCTURAL_ANNOT as S04_STRUCTURAL_ANNOT } from './subworkflows/04_structural_annot'
include { STEP_05_PROTEIN_ALIGNMENT as S05_PROTEIN_ALIGNMENT } from './subworkflows/05_protein_alignment'
include { STEP_06_FUNC_ANNOT as S06_FUNC_ANNOT } from './subworkflows/06_functionnal_annot'
include { STEP_07_TAXO_AFFI as S07_TAXO_AFFI } from './subworkflows/07_taxonomic_affi'
include { STEP_08_BINNING as S08_BINNING } from './subworkflows/08_binning'

include { GET_SOFTWARE_VERSIONS } from './modules/get_software_versions'
include { MULTIQC } from './modules/multiqc'



/*
 * Define helpMessage
 */

 def helpMessage() {

   log.info"""
   Usage:

     The typical command for running the pipeline is as follows:
       nextflow run -profile standard main.nf --input 'samplesheet.csv' --skip_host_filter --skip_kaiju 

     Mandatory arguments:
       --input [path]                Sample sheet: csv file with samples: sample,fastq_1,fastq_2,fasta[for HIFI]
       --type                        Indicate the type of the sequencing data, "SR" : short-read Illumina data use, "HIFI" : long-read PacBio HiFi data 

     Options:
     
     S01_CLEAN_QC options:
       --stop_at_clean               Stop the pipeline at this step
       --skip_clean                  Skip this step.
       --adapter1                    Sequence of adapter 1. Default: Illumina TruSeq adapter.
       --adapter2                    Sequence of adapter 2. Default: Illumina TruSeq adapter.
       --use_sickle                  Use sickle process for high quality trimming.
       --quality_type                Type of quality values for sickle (solexa (CASAVA < 1.3), illumina (CASAVA 1.3 to 1.7), sanger (which is CASAVA >= 1.8)). Default: 'sanger'.
       --skip_host_filter            Skip filter host reads process.
       --host_fasta                  Full path to fasta of host genome ("PATH/name_genome.fasta").
       --host_index                  Full path to directory containing BWA-MEM2 index including base name i.e ("PATH/name_genome.{0123,amb,ann,bwt.2bit.64,pac,sa}").
       You need to use --skip_host_filter or --host_fasta or --skip_01_clean_qc. If it is not the case, an error message will occur.
       --skip_kaiju                  Skip taxonomic affiliation of reads with kaiju.
       --kaiju_verbose               Allow the generation of kaiju verbose output (file can be large)
       --kaiju_db_dir                Directory with kaiju database already built ("PATH/directory").
       --kaiju_db_url                Indicate kaiju database you want to build. Default: "https://kaiju.binf.ku.dk/database/kaiju_db_refseq_2021-02-26.tgz".
       You need to use --kaiju_db_url or --kaiju_db_dir or --skip_kaiju. If it is not the case, an error message will occur.
     
     S02_ASSEMBLY options:
       --stop_at_assembly            Stop the pipeline at this step
       --assembly                    Indicate the assembly tool for short reads ["metaspades" or "megahit" ]. Default: "metaspades".
                                                             or for HiFi reads ["hifiasm-meta", "metaflye"]. Default: "hifiasm-meta".
       --coassembly                  Assemble together samples labeled with the same group in the samplesheet.
              
     S03_FILTERING options:
       --stop_at_filtering           Stop the pipeline at this step
       --skip_filtering              Skip this step
       --min_contigs_cpm [cutoff]    CPM cutoff (Count Per Million) to filter contigs with low number of reads. Default: 1.

     S04_STRUCTURAL_ANNOT options:
       --stop_at_structural_annot    Stop the pipeline at this step

     S05_PROTEIN_ALIGNMENT options:
       --diamond_bank                Path to diamond bank used to align protein sequence of genes: "PATH/bank.dmnd".
                                     This bank must be previously built with diamond makedb.
       
     S06_FUNC_ANNOT options:
       --skip_func_annot             Skip this step
       --percentage_identity [nb]    Sequence identity threshold. Default: 0.95 corresponding to 95%. Use a number between 0 and 1.
       --eggnog_mapper_db_download   Flag --eggnog_mapper_db_download to build the database. Default false: metagWGS didn't build this database.
       --eggnog_mapper_db_dir        Path to eggnog-mapper database "PATH/database_directory/" if it is already built. Default: false.
       You need to use --eggnog_mapper_db_download or --eggnog_mapper_db_dir. If it is not the case, an error message will occur.
       
     S07_TAXO_AFFI options:
       --skip_taxo_affi              Skip this step
       --accession2taxid             Path or FTP adress of file prot.accession2taxid.FULL.gz. Default: "ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz".
       --taxdump                     Path or FTP adress of file taxdump.tar.gz. Default: "ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz".

     SO8_BINNING options:
       --skip_binning                Skip this step
       --gtdbtk_bank                 Path to the GTDBTK database
       --checkm2_bank                Path to the CheckM2 database
       --metabat2_seed               Set the seed for metabat2, for exact reproducibility of metabat2 (default: 0 (random seed))
       --binning_cross_alignment     Mapping strategy to compute co-abundances for binning: .
                                      'individual': Each sample is aligned against the assembly of the sample in question
                                      'group': The samples labeled with the same group in the samplesheet will be aligned against each assembly within the group of samples.
                                      'all': All the samples will be aligned against all the assembly (WARNING: It could take a long time.) (default: individual).
       --min_completeness [nb]       Minimum % of bins completeness for the bins to be kept after bin_refinement step. Default: 50.
       --drep_threshold [nb]         Average Nucleotide Identity (ANI) threshold used for bins de-replication. Default: 0.95 corresponding to 95%. Use a number between 0 and 1.
                                     Most studies agree that 95% ANI is an appropriate threshold for species-level de-replicationton. If the goal of dereplication is to generate a set 
                                     of genomes that are distinct when mapping short reads, 98% (0.98) ANI is an appropriate threshold.

     Other options:
       --outdir                      The output directory where the results will be saved: "dir_name". Default "results".
       --help                        Show this message and exit.
       --multiqc_config              If you want to change the configuration of multiqc report: "PATH/multiqc.yaml". Default "$baseDir/assets/multiqc_config.yaml".

       -profile                      Configuration profile to use.
                                     Available: singularity (use of Singularity containers), conda (use of conda environments),
                                     test and debug profile.
   """.stripIndent()
 }

// Show help message.

if (params.help){
  helpMessage()
  exit 0
}


def getAndCheckHeader() {
    File file = new File(params.input)
    assert file.exists() : "${params.input} file not found"
    def line="";
    file.withReader { reader ->
        line = reader.readLine()
    }
    def tab = line.split(/,/)
    def list = ['sample','flowcell','group','fastq_1','fastq_2', 'assembly']
    for (i in tab) {
        if (!list.contains(i)) {
            exit 1, 'Error 1 while check samplesheet format please enter sample[,flowcell],fastq_1[,fastq_2][,assembly] with header line' 
        }
    } 
    if (!tab.contains("sample") || !tab.contains("fastq_1")){
        exit 1, 'Error 1 while check samplesheet format please enter at least sample,fastq_1 with header line' 
    }
    def header_size = tab.size()
    file.eachLine { row ->
      if (row.size()!=0){
        str = row.split(/,/)
        if (str.size()!=header_size) { 
          exit 1, 'Error 1 while check samplesheet format, each line must have the same number of columns'
        }
        for (val in str){
          if (val.isEmpty()){
            exit 1, 'Error 1 while check samplesheet format, missing elements '
          }
        }
      }
    }
    return tab
}


def returnFile(it) {
    if (it == null) {
        return null
    } else {
        if (!file(it).exists()) exit 1, "Missing file in CSV file: ${it}, see --help for more information"
    }
    return file(it)
}

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}


workflow {
  assembly_tool = params.assembly

  // Check mandatory parameters

  if (params.type.toUpperCase() == 'SR') {

    if ( assembly_tool == null || assembly_tool == ''){
      assembly_tool = 'metaspades'

    }
    else if (!['metaspades','megahit'].contains(assembly_tool)){
      exit 1, "Invalid short read assembler option: ${assembly_tool}. Valid options for short reads: 'metaspades', 'megahit'"
    }

    if (!["solexa","illumina","sanger"].contains(params.quality_type)){
      exit 1, "Invalid quality_type option: ${params.quality_type}. Valid options:'solexa','illumina','sanger'"
    }
    if (!(params.skip_host_filter) && !(params.host_fasta) && !(params.skip_clean)) {
      exit 1, "You must specify --host_fasta or skip cleaning host step with option --skip_host_filter or skip all clean and qc modules with --skip_clean"
    }
  }
  else if (params.type.toUpperCase() == 'HIFI') {

    if ( assembly_tool == null || assembly_tool == ""){
      assembly_tool = 'hifiasm-meta'
    }
    else if (!['hifiasm-meta','metaflye'].contains(assembly_tool)){
      exit 1, "Invalid long read assembler option: ${assembly_tool}. Valid options for HiFi reads: 'hifiasm-meta', 'metaflye'"
    }
  }

//  if ( !(params.stop_at_clean) && !(params.stop_at_assembly) && !(params.stop_at_filtering) && !(params.stop_at_structural_annot) && !(params.diamond_bank) ) {
//      exit 1, "You must specify --stop_at_structural_annot or specify a diamond bank with --diamond_bank"
//  }

  if ( !(params.stop_at_clean) && !(params.stop_at_assembly) && !(params.stop_at_filtering) && !(params.stop_at_structural_annot) && !(params.skip_binning ) && !(params.gtdbtk_bank || params.checkm2_bank) ) {
      exit 1, "You must specify --skip_binning or specify a GTDB-TK bank with --gtdbtk_bank and a checkm2 bank with --checkm2_bank"
  }

  if ( params.coassembly && params.binning_cross_alignment == 'group'){
    exit 1, "--binning_cross_alignment group must not be use use --coassembly."
  }


Channel.fromFilePairs(params.input)
.map { sample,  reads ->
        meta = [
            sample:sample,
        ]
        [ "sample":sample,
          "flowcell":null,
          "group":"all",
          "fastq_1":reads.first(),
          "fastq_2":reads.last(),
          "assembly":params.assembly_result,
          "type":params.type
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
        meta.type = item.type
        meta.assembly = item.assembly
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
      meta.type = item.type
      meta.assembly = item.assembly
      return [meta,item.assembly]
    }
    .set { ch_assembly }


  has_assembly = params.assembly_result //false // (file(params.input).splitCsv ( header:true, sep:',' ).assembly[0] != null)
  has_flowcell = false // (file(params.input).splitCsv ( header:true, sep:',' ).flowcell[0] != null)
  
  ////////////
  // End check samplesheet  
  ////////////

  // Databases
  ch_host_fasta = Channel.empty()
  ch_host_index = Channel.empty()
  ch_kaiju_db = Channel.empty()
  ch_eggnog_db = Channel.empty()
  ch_taxonomy = Channel.empty()
  ch_diamon_db = Channel.empty()
  ch_gtbdtk_db = Channel.empty()
  ch_gtbdtk_db = Channel.empty()

  DATABASES ()
  ch_host_fasta = DATABASES.out.host_fasta
  ch_host_index = DATABASES.out.host_index
  ch_kaiju_db = DATABASES.out.kaiju_db
  ch_eggnog_db = DATABASES.out.eggnog
  ch_taxonomy = DATABASES.out.taxonomy
  ch_diamon_db = DATABASES.out.diamond
  ch_gtbdtk_db = DATABASES.out.gtdbtk
  ch_checkm2_db = DATABASES.out.checkm2

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
  ch_idxstats = Channel.empty()
  ch_unfilter_assembly_flagstat_report = Channel.empty()
  ch_final_assembly_flagstat_report = Channel.empty()
  ch_assembly_report = Channel.empty()
  ch_filtered_report = Channel.empty()
  ch_annotation_report = Channel.empty()
  ch_bins_abundances_report = Channel.empty()
  ch_bins_stats_report = Channel.empty()
  ch_quast = Channel.empty()
  ch_quast_before_filter_report = Channel.empty()
  ch_quast_after_filter_report = Channel.empty()
  ch_software_versions = Channel.empty()
  ch_software_versions_S01 = Channel.empty()
  ch_software_versions_S02 = Channel.empty()
  ch_software_versions_S04 = Channel.empty()
  ch_software_versions_S05 = Channel.empty()
  ch_software_versions_S06 = Channel.empty()
  ch_software_versions_S07 = Channel.empty()
  ch_software_versions_S08 = Channel.empty()
  ch_software_versions_total = Channel.empty()
  ch_circular = Channel.empty()

  ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)

  if ( params.type.toUpperCase() != "SR" &&  params.type.toUpperCase() != "HIFI" )  {
		exit 1, "Invalid type option: ${params.type}. Valid options are 'HiFi' for long-read, 'SR' for short-read"
	}

  if ( !params.skip_clean ) {
    S01_CLEAN_QC (
      ch_reads,
      ch_host_fasta,
      ch_host_index,
      ch_kaiju_db
    )
    ch_reads = S01_CLEAN_QC.out.preprocessed_reads
    ch_cutadapt_report = S01_CLEAN_QC.out.cutadapt_report
    ch_sickle_report = S01_CLEAN_QC.out.sickle_report

    ch_before_host_filter_report = S01_CLEAN_QC.out.before_filter_report

    ch_fastqc_raw_report = S01_CLEAN_QC.out.fastqc_raw_report
    ch_fastqc_clean_report = S01_CLEAN_QC.out.fastqc_clean_report
    ch_kaiju_report = S01_CLEAN_QC.out.kaiju_report

    ch_software_versions_S01 = S01_CLEAN_QC.out.software_versions
  }

  if ( !params.stop_at_clean ) {

    S02_ASSEMBLY ( ch_reads, ch_assembly, has_assembly, assembly_tool, has_flowcell )
    ch_assembly = S02_ASSEMBLY.out.assembly
    ch_reads = S02_ASSEMBLY.out.reads
    ch_bam = S02_ASSEMBLY.out.bam

    ch_sam_coverage = S02_ASSEMBLY.out.coverage
    ch_idxstats = S02_ASSEMBLY.out.idxstats
    ch_unfilter_assembly_flagstat_report = S02_ASSEMBLY.out.flagstat

    ch_quast_before_filter_report = S02_ASSEMBLY.out.assembly_report
    ch_quast = ch_quast_before_filter_report

    ch_software_versions_S02 = S02_ASSEMBLY.out.software_versions
    ch_circular = S02_ASSEMBLY.out.circular
  }

  if ( !params.stop_at_clean && !params.stop_at_assembly && !params.skip_filtering ) {

    ch_min_contigs_cpm = Channel.value(params.min_contigs_cpm)



    S03_FILTERING (
        ch_assembly,
        ch_reads,
        ch_idxstats,
        ch_bam,
        ch_min_contigs_cpm,

    )
    
    ch_assembly = S03_FILTERING.out.selected_contigs
    ch_bam =  S03_FILTERING.out.bam

    

    ch_quast_after_filter_report = S03_FILTERING.out.quast_report
    ch_quast = ch_quast_after_filter_report


    ch_sam_coverage = S03_FILTERING.out.sam_coverage
    ch_final_assembly_flagstat_report = S03_FILTERING.out.sam_flagstat

  }

  ch_annotation_ffn = Channel.empty()
  ch_annotation_faa = Channel.empty()
  ch_annotation_gff = Channel.empty()

  if ( !params.stop_at_clean && !params.stop_at_assembly && !params.stop_at_filtering ) {
    S04_STRUCTURAL_ANNOT ( ch_assembly )
    ch_annotation_ffn = S04_STRUCTURAL_ANNOT.out.ffn
    ch_annotation_faa = S04_STRUCTURAL_ANNOT.out.faa
    ch_annotation_gff = S04_STRUCTURAL_ANNOT.out.gff
    ch_annotation_report = S04_STRUCTURAL_ANNOT.out.report

    ch_software_versions_S04 = S04_STRUCTURAL_ANNOT.out.software_versions
        }

  ch_diamond_result = Channel.empty()

  if ( !params.stop_at_clean && !params.stop_at_assembly && !params.stop_at_filtering && !params.stop_at_structural_annot && (!params.skip_func_annot || !params.skip_taxo_affi)) {
    S05_PROTEIN_ALIGNMENT (ch_annotation_faa, ch_diamon_db)

    ch_diamond_result = S05_PROTEIN_ALIGNMENT.out.diamond_result

    ch_software_versions_S05 = S05_PROTEIN_ALIGNMENT.out.software_versions
  }

  ch_quant_report = Channel.empty()
  ch_v_eggnogmapper = Channel.empty()
  if ( !params.stop_at_clean && !params.stop_at_assembly && !params.stop_at_filtering && !params.stop_at_structural_annot && !params.skip_func_annot) {
      S06_FUNC_ANNOT ( ch_annotation_ffn, ch_annotation_faa, ch_annotation_gff, ch_bam, ch_diamond_result, ch_eggnog_db )
      ch_quant_report = S06_FUNC_ANNOT.out.quant_report

      ch_software_versions_S06 = S06_FUNC_ANNOT.out.software_versions
  }

        if ( !params.stop_at_clean && !params.stop_at_assembly && !params.stop_at_filtering && !params.stop_at_structural_annot && !params.skip_taxo_affi) {
                S07_TAXO_AFFI ( ch_taxonomy, ch_diamond_result, ch_sam_coverage)

    ch_software_versions_S07 = S07_TAXO_AFFI.out.software_versions
        }


  if ( !params.stop_at_clean && !params.stop_at_assembly && !params.stop_at_filtering && !params.stop_at_structural_annot && !params.skip_binning ) {

    S08_BINNING( ch_reads, ch_assembly, ch_bam, ch_gtbdtk_db, ch_checkm2_db, ch_quast, ch_circular)
    ch_bins_abundances_report = S08_BINNING.out.bins_abundances_report
    ch_bins_stats_report = S08_BINNING.out.bins_stats_report

    ch_software_versions_S08 = S08_BINNING.out.software_versions
  }

  ch_software_versions_total = ch_software_versions_S01.mix(ch_software_versions_S02,
                              ch_software_versions_S04,
                              ch_software_versions_S05,
                              ch_software_versions_S06,
                              ch_software_versions_S07,
                              ch_software_versions_S08).unique { it.getBaseName() }.collect()


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
    ch_unfilter_assembly_flagstat_report.collect().ifEmpty([]),
    ch_quast_before_filter_report.collect().ifEmpty([]),
    ch_quast_after_filter_report.collect().ifEmpty([]),
    ch_final_assembly_flagstat_report.collect().ifEmpty([]),
    ch_annotation_report.collect().ifEmpty([]),
    ch_quant_report.collect().ifEmpty([]),
    ch_bins_abundances_report.collect().ifEmpty([]),
    ch_bins_stats_report.collect().ifEmpty([])

  )
  multiqc_report = MULTIQC.out.report

}
