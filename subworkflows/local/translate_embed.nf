include { AMULETY_TRANSLATE  } from '../../modules/nf-core/amulety/translate/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ANTIBERTY} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ANTIBERTA2} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ESM2} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_BALMPAIRED} from '../../modules/nf-core/amulety/embed/main'

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
        AMULETY_EMBED_ANTIBERTY(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain,
            "antiberty"
        )
        ch_versions = ch_versions.mix(AMULETY_EMBED_ANTIBERTY.out.versions)
    }

    if (params.embeddings && params.embeddings.split(',').contains('antiberta2') ){
        AMULETY_EMBED_ANTIBERTA2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain,
            "antiberta2"
        )
        ch_versions = ch_versions.mix(AMULETY_EMBED_ANTIBERTA2.out.versions)
    }

    if (params.embeddings && params.embeddings.split(',').contains('esm2') ){
        AMULETY_EMBED_ESM2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain,
            "esm2"
        )
        ch_versions = ch_versions.mix(AMULETY_EMBED_ESM2.out.versions)
    }

    if (params.embeddings && params.embeddings.split(',').contains('balmpaired') ){
        AMULETY_EMBED_BALMPAIRED(
            AMULETY_TRANSLATE.out.repertoire_translated,
            params.embedding_chain,
            "balm-paired"
        )
        ch_versions = ch_versions.mix(AMULETY_EMBED_BALMPAIRED.out.versions)
    }


    emit:
    versions = ch_versions
}
