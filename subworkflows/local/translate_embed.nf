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
    ch_versions = ch_versions.mix(AMULETY_TRANSLATE.out.versions)

    if (params.embeddings && params.embeddings.split(',').contains('antiberty') ){
        AMULETY_ANTIBERTY(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain
        )
        ch_versions = ch_versions.mix(AMULETY_ANTIBERTY.out.versions)
    }

    if (params.embeddings && params.embeddings.split(',').contains('antiberta2') ){
        AMULETY_ANTIBERTA2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain
        )
        ch_versions = ch_versions.mix(AMULETY_ANTIBERTA2.out.versions)
    }

    if (params.embeddings && params.embeddings.split(',').contains('esm2') ){
        AMULETY_ESM2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain
        )
        ch_versions = ch_versions.mix(AMULETY_ESM2.out.versions)
    }

    if (params.embeddings && params.embeddings.split(',').contains('balmpaired') ){
        AMULETY_BALMPAIRED(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain
        )
        ch_versions = ch_versions.mix(AMULETY_BALMPAIRED.out.versions)
    }


    emit:
    versions = ch_versions
}
