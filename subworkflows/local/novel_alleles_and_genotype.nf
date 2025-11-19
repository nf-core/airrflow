include { NOVEL_ALLELE_INFERENCE } from '../../modules/local/enchantr/novel_allele_inference'
include { BAYESIAN_GENOTYPE_INFERENCE  } from '../../modules/local/enchantr/bayesian_genotype_inference'
include { REASSIGN_ALLELES as REASSIGN_ALLELES_NOVEL; REASSIGN_ALLELES as REASSIGN_ALLELES_GENOTYPE} from '../../modules/local/enchantr/reassign_alleles'


workflow CLONAL_ANALYSIS {
    take:
    ch_repertoire
    ch_reference_fasta
    ch_logo

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()


    // infer novel alleles
    NOVEL_ALLELE_INFERENCE (
        ch_repertoire,
        ch_reference_fasta
    )

    // reassign novel alleles (we can skip this step if no novel alleles were inferred)

    REASSIGN_ALLELES_NOVEL (
        ch_repertoire,
        NOVEL_ALLELE_INFERENCE.out.reference
    )


    // infer genotype (gets the reference from novel alleles inference in any case)

    BAYESIAN_GENOTYPE_INFERENCE (
        REASSIGN_ALLELES_NOVEL.out.repertoire,
        NOVEL_ALLELE_INFERENCE.out.reference
    )

    // reassign genotypes (gets the reference from genotype inference in any case)

    REASSIGN_ALLELES_GENOTYPE (
        REASSIGN_ALLELES_NOVEL.out.repertoire,
        BAYESIAN_GENOTYPE_INFERENCE.out.reference
    )


    emit:
    repertoire = ch_repertoire
    versions = ch_versions
    logs = ch_logs
}

// Function to map
def get_meta_tabs(arr) {
    def meta = [:]
    meta.id            = [arr[0]].unique().join("")
    meta.sample_ids         = arr[1]
    meta.subject_id         = arr[2]
    meta.species            = arr[3]
    meta.single_cell        = arr[4].unique().join("")
    meta.locus              = arr[5].unique().join("")

    def array = []

        array = [ meta, arr[6].flatten() ]

    return array
}
