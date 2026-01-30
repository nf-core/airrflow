include { NOVEL_ALLELE_INFERENCE } from '../../modules/local/enchantr/novel_allele_inference'
include { BAYESIAN_GENOTYPE_INFERENCE  } from '../../modules/local/enchantr/bayesian_genotype_inference'
include { REASSIGN_ALLELES as REASSIGN_ALLELES_NOVEL; REASSIGN_ALLELES as REASSIGN_ALLELES_GENOTYPE} from '../../modules/local/enchantr/reassign_alleles'
include { CLONAL_ANALYSIS } from './clonal_analysis.nf'
include { CLONAL_ASSIGNMENT as CLONAL_ASSIGNMENT_COMPUTE  } from '../../modules/local/enchantr/clonal_assignment'

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

    // infer novel alleles
    if (params.novel_allele_inference) {
        NOVEL_ALLELE_INFERENCE (
            ch_grouped_repertoires,
            ch_reference_fasta
        )

        // reassign novel alleles (we can skip this step if no novel alleles were inferred)

        REASSIGN_ALLELES_NOVEL (
            ch_grouped_repertoires,
            NOVEL_ALLELE_INFERENCE.out.reference,
            ["v"]
        )
        ch_for_genotyping = REASSIGN_ALLELES_NOVEL.out.tab
        ch_for_reference = NOVEL_ALLELE_INFERENCE.out.reference
    } else {
        ch_for_genotyping = ch_grouped_repertoires
        ch_for_reference = ch_reference_fasta
    }

    if (params.single_clone_representative) {
        // TODO: Check if we need the cloneby parameter, or here it can be the same as genotypeby.
        CLONAL_ASSIGNMENT_COMPUTE(
            ch_for_genotyping,
            [params.genotype_clone_threshold],
            ch_reference_fasta.collect(),
            []
        )

        // CLONAL_ANALYSIS(
        //             ch_for_genotyping,
        //             ch_for_reference,
        //             ch_logo.collect().ifEmpty([])
        //         )
        // ch_versions = ch_versions.mix( CLONAL_ANALYSIS.out.versions)

        ch_for_genotyping = CLONAL_ASSIGNMENT_COMPUTE.out.tab//CLONAL_ANALYSIS.out.repertoire
    }

    // infer genotype
    BAYESIAN_GENOTYPE_INFERENCE (
        ch_for_genotyping,
        ch_for_reference
    )

    // reassign genotypes
    REASSIGN_ALLELES_GENOTYPE (
        ch_for_genotyping,
        BAYESIAN_GENOTYPE_INFERENCE.out.reference,
        ["auto"]
    )

    emit:
    repertoire = REASSIGN_ALLELES_GENOTYPE.out.tab
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
