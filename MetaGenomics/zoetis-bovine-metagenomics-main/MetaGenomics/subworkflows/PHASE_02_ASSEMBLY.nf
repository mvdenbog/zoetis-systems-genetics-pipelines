/*--------------------------------------------------

  Author: Mathias Vandenbogaert, December 2023
---------------------------------------------------*/


include { SPADES } from '../modules/assembly'
include { SPADES_KAIJU } from '../modules/assembly'

workflow ASSEMBLY {
	take:
		fq_species	
		kaiju_db
	main:

		SPADES(fq_species)
		ch_assembler_v = SPADES.out.v_spades
		SPADES_KAIJU(SPADES.out.spades_fa_result, kaiju_db)
		ch_kaiju_v = SPADES_KAIJU.out.v_kaiju

		ch_software_versions = ch_assembler_v.first().mix(ch_kaiju_v.first())

	emit:
		software_versions = ch_software_versions
}
