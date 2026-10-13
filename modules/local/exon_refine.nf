process EXON_REFINE {
    tag "${gtf.simpleName}"
    label 'process_low'

    publishDir "${params.outdir}/refined_gtf", mode: params.publish_dir_mode ?: 'copy'

    input:
    path gtf

    output:
    tuple path(gtf), path("${gtf.baseName}_refined.gtf"), emit: refined_gtf

    script:
    def prefix = "${gtf.baseName}_refined"
    """
    OS=\$(uname -s)
    if [ "\$OS" = "Darwin" ]; then
        BINARY="${projectDir}/bin/darwin_x86_64/exonRefine"
    else
        BINARY="${projectDir}/bin/linux_x86_64/exonRefine"
    fi
    \$BINARY ${gtf} -p ${prefix}
    """

    stub:
    """
    touch ${gtf.baseName}_refined.gtf
    """
}
