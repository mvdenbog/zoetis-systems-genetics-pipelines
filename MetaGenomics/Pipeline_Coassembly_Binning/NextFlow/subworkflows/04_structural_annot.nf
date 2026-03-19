include { PRODIGAL } from '../modules/prodigal'
include { BARRNAP } from '../modules/barrnap'
include { TRNASCAN_SE } from '../modules/trnascan_se'
include { MERGE_ANNOTATIONS } from '../modules/merge_annotations'


workflow STEP_04_STRUCTURAL_ANNOT {
  take: assembly

  main:
    PRODIGAL( assembly )
    BARRNAP( assembly )
    TRNASCAN_SE( assembly )

    annotations_ch = assembly.join(PRODIGAL.out.faa).join(PRODIGAL.out.gff).join(BARRNAP.out.gff)
                                              .join(TRNASCAN_SE.out.gff)

    MERGE_ANNOTATIONS(annotations_ch)

    ch_software_versions = PRODIGAL.out.v_prodigal.first().mix( BARRNAP.out.v_barrnap.first(), 
                                                                TRNASCAN_SE.out.v_tRNAscan.first())

  emit:
    report = MERGE_ANNOTATIONS.out.report
    ffn = MERGE_ANNOTATIONS.out.ffn
    gff = MERGE_ANNOTATIONS.out.gff
    faa = MERGE_ANNOTATIONS.out.faa
    software_versions = ch_software_versions
}