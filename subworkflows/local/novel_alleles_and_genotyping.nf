include { NOVEL_ALLELE_INFERENCE } from '../../modules/local/enchantr/novel_allele_inference'
include { BAYESIAN_GENOTYPE_INFERENCE  } from '../../modules/local/enchantr/bayesian_genotype_inference'
include { REASSIGN_ALLELES as REASSIGN_ALLELES_NOVEL; REASSIGN_ALLELES as REASSIGN_ALLELES_GENOTYPE} from '../../modules/local/enchantr/reassign_alleles'
include { CLONAL_ANALYSIS } from './clonal_analysis.nf'

workflow NOVEL_ALLELES_AND_GENOTYPING {
    take:
    ch_repertoire
    ch_reference_fasta
    ch_validated_samplesheet
    ch_logo

    main:
    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    // merge all repertoires by genotypeby metadata field
    ch_repertoire.map{ it -> [ it[0]."${params.genotypeby}",
                                it[0].id,
                                it[0].subject_id,
                                it[0].species,
                                it[0].single_cell,
                                it[0].locus,
                                it[1] ] }
                .groupTuple()
                .map{ get_meta_tabs(it) }
                .set{ ch_grouped_repertoires }

    //TODO: conditional on params.novel_allele_inference
    // infer novel alleles
    if (params.novel_allele_inference) {
        NOVEL_ALLELE_INFERENCE (
            ch_grouped_repertoires,
            ch_reference_fasta,
            ch_validated_samplesheet.collect()
        )

        // reassign novel alleles (we can skip this step if no novel alleles were inferred)

        REASSIGN_ALLELES_NOVEL (
            ch_grouped_repertoires,
            NOVEL_ALLELE_INFERENCE.out.reference,
            ch_validated_samplesheet.collect(),
            "v" //TODO: update this to pass actual segments. We only need to reassign V after novel allele inference.
        )
    }

    // TODO: what are we doing with the reference if we are not running novel allele inference?

    // infer clones (gets the reference from novel alleles inference in any case)

    CLONAL_ANALYSIS(
                REASSIGN_ALLELES_NOVEL.out.repertoire,
                NOVEL_ALLELE_INFERENCE.out.reference,
                ch_logo.collect().ifEmpty([])
            )
    ch_versions = ch_versions.mix( CLONAL_ANALYSIS.out.versions)

    // infer genotype (gets the reference from novel alleles inference in any case)

    BAYESIAN_GENOTYPE_INFERENCE (
        REASSIGN_ALLELES_NOVEL.out.repertoires,
        NOVEL_ALLELE_INFERENCE.out.reference,
        ch_validated_samplesheet.collect()
    )

    // reassign genotypes (gets the reference from genotype inference in any case)

    REASSIGN_ALLELES_GENOTYPE (
        REASSIGN_ALLELES_NOVEL.out.repertoires,
        BAYESIAN_GENOTYPE_INFERENCE.out.reference,
        ch_validated_samplesheet.collect(),
        "auto" //TODO: update this to pass actual segments. We're running over all segment after genotype inference.
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
