// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MERGE_UMI {
    tag "$meta.id"

    input:
    tuple val(meta), path(R1), path(R2), path(I1)

    output:
    tuple val(meta), path($R1), path($R2)   , emit: reads

    script:
    if (params.index_file) {
        """
        merge_R1_umi.py -R1 "${R1}" -I1 "${I1}" -o UMI_R1.fastq.gz --umi_start $params.umi_start --umi_length $params.umi_length
        gunzip -f "UMI_R1.fastq.gz" 
        mv "UMI_R1.fastq" "${meta.id}_${R1.baseName}"
        gunzip -f "${R2}"
        mv "${R2.baseName}" "${meta.id}_${R2.baseName}"
        """
    } else {
        """
        gunzip -f "${R1}"
        mv "${R1.baseName}" "${meta.id}_${R1.baseName}"
        gunzip -f "${R2}"
        mv "${R2.baseName}" "${meta.id}_${R2.baseName}"
        """
    }
}