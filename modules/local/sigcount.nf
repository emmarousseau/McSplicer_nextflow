process SIGCOUNT {
    tag "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir}/signature_counts", mode: params.publish_dir_mode ?: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(refined_gtf)

    output:
    tuple val(meta), path("${meta.id}.cnt"), path(refined_gtf), emit: counts

    script:
    def binary_name = meta.paired_end ? 'sigcount_pe' : 'sigcount'
    """
    OS=\$(uname -s)
    if [ "\$OS" = "Darwin" ]; then
        BINARY="${projectDir}/bin/darwin_x86_64/${binary_name}"
    else
        BINARY="${projectDir}/bin/linux_x86_64/${binary_name}"
    fi
    \$BINARY ${bam} ${refined_gtf} ${meta.id}
    """

    stub:
    """
    touch ${meta.id}.cnt
    """
}
