include { DIAMOND } from '../modules/diamond'

workflow STEP_05_PROTEIN_ALIGNMENT {
    take:
        prokka_faa
        diamond

	main:
        ch_diamond_v = Channel.empty()


        ch_diamond_result =Channel.empty()
        
        DIAMOND (
            prokka_faa,
            diamond
        )

        ch_diamond_result = DIAMOND.out.m8
        ch_diamond_v = DIAMOND.out.v_diamond
        
    emit:
        diamond_result = ch_diamond_result
        software_versions = ch_diamond_v.first()
    }
