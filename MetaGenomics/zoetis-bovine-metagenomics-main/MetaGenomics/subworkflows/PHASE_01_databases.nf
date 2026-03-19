
/*--------------------------------------------------
  Author: Mathias Vandenbogaert, December 2023
---------------------------------------------------*/


include { GET_DB_VERSIONS } from '../modules/get_db_versions.nf'

workflow DATABASES {

    main:

        ch_host_fasta = Channel.empty()
        ch_host_index = Channel.empty()
        ch_kaiju_db = Channel.empty()
        ch_kraken2_db = Channel.empty()
        ch_k2_taxonomy_db = Channel.empty()

        if ( !params.skip_clean && !params.skip_host_filter && (params.host_fasta == "" || !params.host_fasta)) {
            exit 1, "You must specify a host fasta file with --host_fasta or skip the host filtering step with --skip_host_filter or --skip_clean"
        }
        
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

        ch_kraken2_db = Channel.value(params.kraken2_db_dir)
        ch_k2_taxonomy_db = Channel.value(params.k2_taxonomy_db_dir)

        ch_taxonomy = Channel.empty()
        ch_taxdump = Channel.empty()
        ch_accession2taxid = Channel.empty()

        GET_DB_VERSIONS(
            ch_host_fasta.ifEmpty([]),
            ch_kaiju_db.map {it -> it[0]}.ifEmpty([]),
            ch_kraken2_db.ifEmpty([]),
            ch_accession2taxid.ifEmpty([]),
            ch_taxdump.ifEmpty([]),
        )
        
    emit:
        host_fasta = ch_host_fasta
        host_index = ch_host_index
        kaiju_db = ch_kaiju_db
        kraken2_db = ch_kraken2_db
        k2_taxonomy_db = ch_k2_taxonomy_db
        taxonomy = ch_taxonomy.first()
}

process INDEX_HOST {
disk 150.GB
memory 100.GB
label "GENERAL"
  publishDir "${params.databases}/IndexHost"

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
  publishDir "${params.databases}/KaijuDB"

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


