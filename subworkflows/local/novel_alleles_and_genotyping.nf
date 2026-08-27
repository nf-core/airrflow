include { NOVEL_ALLELE_INFERENCE } from '../../modules/local/enchantr/novel_allele_inference'
include { BAYESIAN_GENOTYPE_INFERENCE as BAYESIAN_GENOTYPE_INFERENCE_IG } from '../../modules/local/enchantr/bayesian_genotype_inference'
include { BAYESIAN_GENOTYPE_INFERENCE as BAYESIAN_GENOTYPE_INFERENCE_TR } from '../../modules/local/enchantr/bayesian_genotype_inference'
include { REASSIGN_ALLELES as REASSIGN_ALLELES_NOVEL; REASSIGN_ALLELES as REASSIGN_ALLELES_GENOTYPE} from '../../modules/local/enchantr/reassign_alleles'
include { CLONAL_ANALYSIS } from './clonal_analysis'
include { CLONAL_ASSIGNMENT as CLONAL_ASSIGNMENT_GENOTYPING } from '../../modules/local/enchantr/clonal_assignment'

workflow NOVEL_ALLELES_AND_GENOTYPING {
    take:
    ch_repertoire
    ch_reference_fasta
    ch_validated_samplesheet
    ch_logo
    genotypeby
    novel_allele_inference
    single_clone_representative
    genotyping_clonal_threshold
    cloneby
    singlecell

    main:
    ch_logs = channel.empty()

    // Flatten each repertoire into a tuple keyed by the genotypeby field and locus.
    ch_repertoire
        .combine(ch_reference_fasta)
        .map{ it ->
                def meta = it[0]
                def rep = it[1]
                def ref = it[2]
                def genotypeby_field = genotypeby=="sample_id" ? "id" : genotypeby
                def genotype_id = meta[genotypeby_field]
                def locus = meta.locus.toUpperCase()
                [ genotype_id,
                                    locus,
                                    meta.id,
                                    meta.sample_id,
                                    meta.subject_id,
                                    meta.species,
                                    meta.single_cell,
                                    meta.locus,
                                    rep,
                                    ref ] }
                    .set{ ch_repertoires_for_grouping }

    ch_repertoires_for_grouping
                    .groupTuple(by: [0,1])
                    .map{ get_meta_tabs(it) }
                    .set{ ch_grouped_repertoires }

    ch_grouped_repertoires
        .branch { it ->
            ig: it[0].locus.toUpperCase() == 'IG'
            tr: it[0].locus.toUpperCase() == 'TR'
        }
        .set{ ch_grouped_repertoires_by_locus }

    // infer novel alleles for IG only
    if (novel_allele_inference) {
        NOVEL_ALLELE_INFERENCE (
            ch_grouped_repertoires_by_locus.ig
        )

        // reassign novel alleles (we can skip this step if no novel alleles were inferred)
        ch_grouped_repertoires_by_locus.ig
            .join(NOVEL_ALLELE_INFERENCE.out.reference)
            .map { it ->
                def meta = it[0]
                def reps = it[1]
                def new_ref = it[3]
                [ meta, reps, new_ref ]
            }
            .set{ ch_reassign_alleles }

        REASSIGN_ALLELES_NOVEL (
            ch_reassign_alleles,
            ["v"],
            genotypeby //TODO: @ayeletperes check if this is correct
        )

        REASSIGN_ALLELES_NOVEL.out.tab.dump(tag: "reassign alleles novel")

        REASSIGN_ALLELES_NOVEL.out.tab
            .join(NOVEL_ALLELE_INFERENCE.out.reference)
            .mix(ch_grouped_repertoires_by_locus.tr)
            .set{ ch_repertoire_reference }

    } else {
        ch_repertoire_reference = ch_grouped_repertoires
    }

    ch_repertoire_reference
        .branch { it ->
            ig: it[0].locus.toUpperCase() == 'IG'
            tr: it[0].locus.toUpperCase() == 'TR'
        }
        .set{ ch_repertoire_reference_by_locus }

    if (single_clone_representative) {
        // create separate channels for repertoire and reference based on the genotypeby metadata field

        CLONAL_ASSIGNMENT_GENOTYPING(
            ch_repertoire_reference_by_locus.ig,
            [genotyping_clonal_threshold],
            [],
            cloneby,
            singlecell
        )
        CLONAL_ASSIGNMENT_GENOTYPING.out.tab
            .join(ch_repertoire_reference_by_locus.ig
                        .map{ it -> [it[0], it[2]] })
            .set{ ch_ig_for_genotyping }
    } else {
        ch_ig_for_genotyping = ch_repertoire_reference_by_locus.ig
    }

    // infer genotype. TR skips single-clone representative preprocessing.
    BAYESIAN_GENOTYPE_INFERENCE_IG (
        ch_ig_for_genotyping,
        genotypeby,
        single_clone_representative
    )
    BAYESIAN_GENOTYPE_INFERENCE_TR (
        ch_repertoire_reference_by_locus.tr,
        genotypeby,
        false
    )

    // Combine the personalized genotype references inferred for IG and TR loci.
    BAYESIAN_GENOTYPE_INFERENCE_IG.out.reference
        .mix(BAYESIAN_GENOTYPE_INFERENCE_TR.out.reference)
        .set{ ch_genotype_reference }

    // Reassign genotype calls using the personalized genotype reference for both
    // IG and TR loci. The full (pre-single-clone) repertoire is reassigned.
    ch_repertoire_reference_by_locus.ig
        .mix(ch_repertoire_reference_by_locus.tr)
        .map{ it -> [it[0], it[1]] }
        .join(ch_genotype_reference)
        .set{ ch_for_reassign }

    REASSIGN_ALLELES_GENOTYPE (
        ch_for_reassign,
        ["auto"],
        genotypeby
    )

    REASSIGN_ALLELES_GENOTYPE.out.tab.dump(tag: "reassign alleles genotype out tab")

    ch_repertoire_reference = REASSIGN_ALLELES_GENOTYPE.out.tab.join(ch_genotype_reference)
    ch_repertoire_reference.dump(tag: "ch_repertoire_reference_genotyping")


    emit:
    repertoire_reference = ch_repertoire_reference
    logs = ch_logs
}

// Function to map
def get_meta_tabs(arr) {
    def genotype_id = arr[0]
    def grouping_locus = arr[1]

    if (!['IG', 'TR'].contains(grouping_locus)) {
        error "Unsupported locus '${grouping_locus}' found for ${genotype_id}. Genotyping supports IG and TR loci only."
    }

    if (arr[4].unique().size() > 1) {
        error "Multiple subject IDs found for ${genotype_id} (${arr[4].join(', ')}). It is not possible to perform joint genotyping of samples from different subjects. Please check the 'genotypeby' parameter."
    }

    def meta = [:]
    meta.id                 = "${genotype_id}_${grouping_locus}"
    meta.sample_id          = arr[3].flatten()
    meta.subject_id         = arr[4].unique().join("")
    meta.species            = arr[5].unique().join("")
    meta.single_cell        = arr[6].unique().join("")
    meta.locus              = grouping_locus
    def array = []

    array = [ meta, arr[8].flatten(), arr[9].unique() ]
    if (arr[9].unique().size() > 1) {
        error "Multiple reference fasta files found for ${meta.id}."
    }


    return array
}
