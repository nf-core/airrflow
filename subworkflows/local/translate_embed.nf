include { AMULETY_TRANSLATE  } from '../../modules/nf-core/amulety/translate/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ANTIBERTY} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ANTIBERTA2} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_ESM2} from '../../modules/nf-core/amulety/embed/main'
include { AMULETY_EMBED  as AMULETY_EMBED_BALMPAIRED} from '../../modules/nf-core/amulety/embed/main'

workflow TRANSLATE_EMBED {
    take:
    ch_repertoire
    ch_reference_igblast
    embeddings
    embedding_chain

    main:

    AMULETY_TRANSLATE(
        ch_repertoire,
        ch_reference_igblast
    )

    if (embeddings && embeddings.split(',').contains('antiberty') ){
        AMULETY_EMBED_ANTIBERTY(
            AMULETY_TRANSLATE.out.repertoire_translated,
            embedding_chain,
            "antiberty"
        )
    }

    if (embeddings && embeddings.split(',').contains('antiberta2') ){
        AMULETY_EMBED_ANTIBERTA2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            embedding_chain,
            "antiberta2"
        )
    }

    if (embeddings && embeddings.split(',').contains('esm2') ){
        AMULETY_EMBED_ESM2(
            AMULETY_TRANSLATE.out.repertoire_translated,
            embedding_chain,
            "esm2"
        )
    }

    if (embeddings && embeddings.split(',').contains('balmpaired') ){
        AMULETY_EMBED_BALMPAIRED(
            AMULETY_TRANSLATE.out.repertoire_translated,
            embedding_chain,
            "balm-paired"
        )
    }

}
