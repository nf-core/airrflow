include { CELLRANGER_MKVDJREF  } from '../../modules/local/cellranger/mkvdjref'
include { CELLRANGER_VDJ    } from '../../modules/nf-core/cellranger/vdj/main'

workflow SC_RAW_INPUT {

    take:
    

    main:

    ch_versions = Channel.empty()
    ch_logs = Channel.empty()

    // validate library generation method parameter
    if (params.library_generation_method == 'specific_pcr_5p_race_umi') {
        if (params.vprimers) {
            error "The transcript-specific primer, 5'-RACE, UMI library generation method does not require V-region primers, please provide a reference file instead or select another library method option."
        } else if (params.race_linker) {
            error "The transcript-specific primer, 5'-RACE, UMI library generation method does not require the --race_linker parameter, please provide a reference file instead or select another library method option."
        } 
        if (params.cprimers)  {
            error "The transcript-specific primer, 5'-RACE, UMI library generation method does not require C-region primers, please provide a reference file instead or select another library method option."
        }
        if (params.umi_length > 0)  {
            error "The transcript-specific primer, 5'-RACE, UMI library generation method does not require to set the UMI length, please provide a reference file instead or select another library method option."
        } 
        if (params.sc_reference)  {
            ch_sc_refence = Channel.fromPath(params.sc_reference, checkIfExists: true)
        } else {
            error "The transcript-specific primer, 5'-RACE, UMI library generation method requires to provide a reference file."
        }
    }


}
