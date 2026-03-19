include { GET_DB_VERSIONS } from '../modules/get_db_versions.nf'

workflow DATABASES {

    main:

        ch_host_fasta = Channel.empty()
        ch_host_index = Channel.empty()
        ch_kaiju_db = Channel.empty()

        if ( !params.skip_clean && !params.skip_host_filter && (params.host_fasta == "" || !params.host_fasta)) {
            exit 1, "You must specify a host fasta file with --host_fasta or skip the host filtering step with --skip_host_filter or --skip_clean"
        }
        
        if ( params.type.toUpperCase() == "SR" ) {
            if ( !params.skip_clean && !params.skip_host_filter ) {
                ch_host_fasta = Channel.value(file(params.host_fasta, checkIfExists: true))
                if ( params.host_index == "" ||!params.host_index ) {
                    INDEX_HOST(ch_host_fasta)
                    ch_host_index = INDEX_HOST.out.index
                }
                else {
                    ch_host_index = Channel.value(file(params.host_index))
                }
            }
        }
        else if ( params.type.toUpperCase() == "HIFI" )  {
            if ( !params.skip_clean && !params.skip_host_filter ) {
                ch_host_fasta = Channel.value(file(params.host_fasta, checkIfExists: true))
            }

        }


        
        if ( !params.skip_clean && !params.skip_kaiju ) { //kaiju_db
            if ( !params.kaiju_db_dir && params.kaiju_db_url ) {
                INDEX_KAIJU(params.kaiju_db_url)
                ch_kaiju_db = INDEX_KAIJU.out.kaiju_db
            } else if (params.kaiju_db_dir) {
                    if (file(params.kaiju_db_dir + "/kaiju_db*.fmi").size == 1) {
                        ch_kaiju_db = Channel.value([file(params.kaiju_db_dir + "/nodes.dmp"), file(params.kaiju_db_dir + "/kaiju_db*.fmi"), file(params.kaiju_db_dir + "/names.dmp")])
                    } 
                    else if (!file(params.kaiju_db_dir).isDirectory()) {
                        exit 1, "kaiju_db_dir ${params.kaiju_db_dir} does not exists."
                    }
                    else if (file(params.kaiju_db_dir + "/kaiju_db*.fmi").size > 1) {
                        exit 1, "There is more than one file ending with .fmi in kaiju_db_dir: ${params.kaiju_db_dir}"
                    }
                    else if (file(params.kaiju_db_dir + "/kaiju_db*.fmi").size == 0) {
                        exit 1, "There is no file ending with .fmi in kaiju_db_dir: ${params.kaiju_db_dir}"
                    }
            } else {
                exit 1, "You must specify --kaiju_db_url or --kaiju_db_dir"
            }
        }

        ch_eggnog = Channel.empty()
        if ( !params.stop_at_clean && !params.stop_at_filtering && !params.stop_at_assembly && !params.stop_at_structural_annot && !params.skip_func_annot ) { //eggnog_mapper_db
            if( params.eggnog_mapper_db_dir != "" ) {
                ch_eggnog = Channel.fromPath(params.eggnog_mapper_db_dir, checkIfExists: true).first()
            }
            else if ( params.eggnog_mapper_db_download ) {
                // Built eggNOG-mapper database.
                EGGNOG_MAPPER_DB()
                ch_eggnog = EGGNOG_MAPPER_DB.out.functional_annot_db
            } else {
                exit 1, "You must specify --eggnog_mapper_db_download or --eggnog_mapper_db_dir"
            }
        }

        ch_taxonomy = Channel.empty()
        ch_taxdump = Channel.empty()
        ch_accession2taxid = Channel.empty()
        if ( !params.stop_at_clean && !params.stop_at_filtering && !params.stop_at_assembly && !params.stop_at_structural_annot && !params.skip_taxo_affi ) {
        
            ch_taxdump = Channel.fromPath(file(params.taxdump), checkIfExists: true)
            
                    
            ch_accession2taxid = Channel.fromPath(file(params.accession2taxid), checkIfExists: true)
            
            ch_taxonomy = ch_accession2taxid.combine(ch_taxdump)
               
        }

        ch_diamond = Channel.empty()
        if ( !(params.stop_at_clean) && !(params.stop_at_assembly) && !(params.stop_at_filtering) && !(params.stop_at_structural_annot) && (!(params.skip_taxo_affi) || !(params.skip_func_annot)) ) {
            ch_diamond =Channel.value(file(params.diamond_bank))
        }

        ch_gtdbtk_db = Channel.empty()
        ch_checkm2_db = Channel.empty()

        if ( !(params.stop_at_clean) && !(params.stop_at_assembly) && !(params.stop_at_filtering) && !(params.stop_at_structural_annot) && !(params.skip_binning) ) {
            ch_gtdbtk_db = Channel.fromPath(params.gtdbtk_bank).first()
            ch_checkm2_db = Channel.fromPath(params.checkm2_bank).first()
        }
        
        GET_DB_VERSIONS(
            ch_host_fasta.ifEmpty([]),
            ch_kaiju_db.map {it -> it[0]}.ifEmpty([]),
            ch_eggnog.ifEmpty([]),
            ch_accession2taxid.ifEmpty([]),
            ch_taxdump.ifEmpty([]),
            ch_diamond.ifEmpty([]),
            ch_gtdbtk_db.ifEmpty([]),
            ch_checkm2_db.ifEmpty([])
        )
        
    emit:
        host_fasta = ch_host_fasta
        host_index = ch_host_index
        kaiju_db = ch_kaiju_db
        eggnog = ch_eggnog
        taxonomy = ch_taxonomy.first()
        diamond = ch_diamond
        gtdbtk = ch_gtdbtk_db
        checkm2 = ch_checkm2_db
}

process INDEX_HOST {
  publishDir "${params.databases}/index_host"

  input:
  path fasta

  output:
  path "${fasta}.*", emit: index

  script:
  """
  bwa-mem2 index $fasta
  """
}

process INDEX_KAIJU {
  publishDir "${params.databases}/kaiju_db"

  input:
    val database

  output:
    tuple path("nodes.dmp"), path("*.fmi"), path("names.dmp"), emit: kaiju_db

  script:
    """
    wget ${database}
    file='${database}'
    fileNameDatabase=\${file##*/}
    echo \$fileNameDatabase
    tar -zxvf \$fileNameDatabase
    """
}

process EGGNOG_MAPPER_DB {
    publishDir "${params.databases}/eggnog_db"
    label 'EGGNOG'

    output:
        path "db_eggnog_mapper" , emit:  functional_annot_db

    script:
        """
        mkdir db_eggnog_mapper
        download_eggnog_data.py -f -y --data_dir db_eggnog_mapper
        """
}

