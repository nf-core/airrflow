include { AMULETY_TRANSLATE } from '../../modules/local/amulety/translate/main.nf'
include { AMULETY_ANTIBERTY } from '../../modules/local/amulety/antiberty/main.nf'

workflow TRANSLATE_EMBED {
    take:
    ch_repertoire
    ch_reference_igblast

    main:
    ch_versions = Channel.empty()

    AMULETY_TRANSLATE(
        ch_repertoire,
        ch_reference_igblast
    )

    if (params.embeddings && params.embeddings.split(',').contains('antiberty') ){
        AMULETY_ANTIBERTY(
            AMULETY_TRANSLATE.out.repertoire_translated
        )
    }

    emit:
    versions = ch_versions
}
