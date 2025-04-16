include { AMULETY_TRANSLATE } from '../../modules/nf-core/amulety/translate/main'
include { AMULETY_ANTIBERTY } from '../../modules/nf-core/amulety/antiberty/main'

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
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain
        )
    }

    emit:
    versions = ch_versions
}
