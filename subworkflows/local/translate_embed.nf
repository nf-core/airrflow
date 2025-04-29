include { AMULETY_TRANSLATE  } from '../../modules/nf-core/amulety/translate/main'
include { AMULETY_ANTIBERTY  } from '../../modules/nf-core/amulety/antiberty/main'
include { AMULETY_ANTIBERTA2 } from '../../modules/nf-core/amulety/antiberta2/main'
include { AMULETY_BALMPAIRED } from '../../modules/nf-core/amulety/balmpaired/main'
include { AMULETY_ESM2       } from '../../modules/nf-core/amulety/esm2/main'

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

    if (params.embeddings && params.embeddings.split(',').contains('antiberta2') ){
        AMULETY_ANTIBERTA2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain
        )
    }

    if (params.embeddings && params.embeddings.split(',').contains('esm2') ){
        AMULETY_ESM2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain
        )
    }

    if (params.embeddings && params.embeddings.split(',').contains('balmpaired') ){
        AMULETY_BALMPAIRED(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain
        )
    }


    emit:
    versions = ch_versions
}
