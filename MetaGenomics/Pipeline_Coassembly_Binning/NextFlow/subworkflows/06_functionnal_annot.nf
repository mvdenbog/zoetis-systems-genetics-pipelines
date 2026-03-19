include { CD_HIT } from '../modules/cd_hit'
include { QUANTIFICATION } from '../modules/feature_counts'
include { EGGNOG_MAPPER } from '../modules/eggnog_mapper'
include { BEST_HITS } from '../modules/best_hits'
include { MERGE_QUANT_ANNOT_BEST } from '../modules/merge_quant_eggnog_best'
include { FUNCTIONAL_ANNOT_TABLE } from '../modules/functional_annot_table'

// cd_hit + quantification + quantification_table + eggnog_mapper_db + eggnog_mapper 
// + best_hit_diamond + merge_quantif and functionnal annot + make_functionnal_annotation_tables

workflow STEP_06_FUNC_ANNOT {
    take:
        ffn // channel: [ val(meta), path(ffn) ]
        faa // channel: [ val(meta), path(faa) ]
        gff // channel: [ val(meta), path(gff) ]
        bam // channel: [ val(meta), path(bam), path(bam_index) ]
        m8 // channel: [ val(meta), path(diamond_file) ]
        eggnog_db

    main:

        ch_cdhit_v = Channel.empty()
        ch_featurecounts_v = Channel.empty()
        ch_eggnogmapper_v = Channel.empty()
    
        CD_HIT ( faa, params.percentage_identity )
        ch_individual_clstr_table = CD_HIT.out.individual_clstr_table
        ch_global_clstr_table = CD_HIT.out.global_clstr_table
        ch_cdhit_v = CD_HIT.out.v_cdhit

        if (params.coassembly){
            bam.map { meta, bam, bai ->
                        [ meta.group, meta, bam, bai ]}
                    .set { ch_bam_tmp }

            gff.map { meta, gff ->
                        [ meta.group, gff]} 
                    .combine( ch_bam_tmp, by: 0)
                    .map { group, gff, meta, bam, bai -> 
                        [ meta, gff, bam, bai ]}
                    .set { ch_gff_and_bam }
        } else {
            ch_gff_and_bam = gff.join(bam, remainder: false)
        }

        QUANTIFICATION ( ch_gff_and_bam, ch_individual_clstr_table, ch_global_clstr_table)
        ch_quant_table = QUANTIFICATION.out.quantification_table
        ch_quant_report = QUANTIFICATION.out.quant_report
        ch_featurecounts_v = QUANTIFICATION.out.v_featurecounts

        EGGNOG_MAPPER ( faa, eggnog_db )
        ch_annot = EGGNOG_MAPPER.out.annot.collect()
        ch_eggnogmapper_v = EGGNOG_MAPPER.out.version

        BEST_HITS ( m8 )
        ch_best_hits = BEST_HITS.out.best_hits.collect()

        MERGE_QUANT_ANNOT_BEST ( ch_quant_table, ch_annot, ch_best_hits )
        ch_merged_quant_annot_best = MERGE_QUANT_ANNOT_BEST.out.merged

        FUNCTIONAL_ANNOT_TABLE ( ch_merged_quant_annot_best )

        ch_software_versions = ch_cdhit_v.first().mix(ch_featurecounts_v.first(),
                                                       ch_eggnogmapper_v.first())
        

    emit:
        functional_annot = FUNCTIONAL_ANNOT_TABLE.out.functional_annot
        quant_report = ch_quant_report
        software_versions = ch_software_versions
}
