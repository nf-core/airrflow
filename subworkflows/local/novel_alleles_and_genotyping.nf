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
    def outputby = params.genotypeby=="sample_id" ? "id" : params.genotypeby //TODO: we need to change this so we can handle the cases of inferring based on naive and reassigning all
    // merge all repertoires by genotypeby metadata field
    ch_repertoire
        .combine(ch_reference_fasta)
        .map{ it -> 
             def meta = it[0]
             def rep = it[1]
             def ref = it[2]
             def genotypeby = params.genotypeby=="sample_id" ? "id" : params.genotypeby
             [ meta."${genotypeby}",
                                meta.id,
                                meta.subject_id,
                                meta.species,
                                meta.single_cell,
                                meta.locus,
                                rep,
                                ref ] }
                .groupTuple()
                .map{ get_meta_tabs(it) }
                .set{ ch_grouped_repertoires }

    // infer novel alleles
    if (params.novel_allele_inference) {
        NOVEL_ALLELE_INFERENCE (
            ch_grouped_repertoires
        )

        // reassign novel alleles (we can skip this step if no novel alleles were inferred)
        ch_grouped_repertoires
            .join(NOVEL_ALLELE_INFERENCE.out.reference)
            .map { it -> 
                def meta = it[0]
                def reps = it[1]
                def new_ref = it[3]
                [ meta, reps, new_ref ]
            }
            .set{ ch_for_genotyping }
        
        REASSIGN_ALLELES_NOVEL (
            ch_for_genotyping,
            ["v"],
            outputby
        )

        REASSIGN_ALLELES_NOVEL.out.tab
            .join(NOVEL_ALLELE_INFERENCE.out.reference)
            .set{ ch_for_genotyping }
        

    } else {
        ch_for_genotyping = ch_grouped_repertoires
    }

    if (params.single_clone_representative) {
        // TODO: Check if we need the cloneby parameter, or here it can be the same as genotypeby.
        // create separate channels for repertoire and reference based on the genotypeby metadata field
        ch_for_genotyping
            .map{ it -> [it[0], it[1]] }
            .set{ ch_for_genotyping_rep }
        ch_for_genotyping
            .map{ it -> it[2] }
            .set{ ch_for_genotyping_ref }
        CLONAL_ASSIGNMENT_COMPUTE(
            ch_for_genotyping_rep,
            [params.genotype_clone_threshold],
            ch_for_genotyping_ref,
            []
        )
        CLONAL_ASSIGNMENT_COMPUTE.out.tab
            .join(ch_for_genotyping
            .map{ it -> [it[0], it[2]] })
            .set{ ch_for_genotyping }
    }

    // infer genotype
    BAYESIAN_GENOTYPE_INFERENCE (
        ch_for_genotyping
    )
    
    ch_grouped_repertoires
        .map{ it -> [it[0], it[1]] }
        .join(BAYESIAN_GENOTYPE_INFERENCE.out.reference)
        .set{ ch_for_reassign }

    // reassign genotypes
    REASSIGN_ALLELES_GENOTYPE (
        ch_for_reassign,
        ["auto"],
        outputby
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

    array = [ meta, arr[6].flatten(), arr[7][0] ]
    return array
}
