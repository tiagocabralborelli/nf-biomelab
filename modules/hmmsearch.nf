process ANNOTAION {
    publishDir 'results/hmmsearch/', mode:'symlink'
    conda '${params.envs.training}'
    input:
        tuple val(sample_name), path(nucl), path(prot)
        val(hmm_args)
        val(tag)

    output:
        tuple val(sample_name), path("${sample_name}_${tag}.tblout")

    script:
    """
    hmmsearch --cpu 12 --tblout ${sample_name}_${tag}.tblout ${hmm_args} ${prot}

    """
}