process MCSPLICER {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outdir}/mcsplicer", mode: params.publish_dir_mode ?: 'copy'

    input:
    tuple val(meta), path(refined_gtf), path(cnt_file)

    output:
    tuple val(meta), path("output/**/*.csv"), emit: results

    script:
    def gene_id_arg   = meta.gene_id    ? "--gene_id ${meta.gene_id}"          : ''
    def gene_list_arg = params.gene_list ? "--gene_list ${params.gene_list}"    : ''
    def prefix_arg    = meta.prefix     ? "--prefix ${meta.prefix}"            : ''
    def bootstrap_arg = meta.bootstraps.toInteger() > 0 ? "--bootstraps ${meta.bootstraps}" : ''
    def anno_arg      = meta.anno == 'y' ? "--anno y"                          : ''
    """
    python3 ${projectDir}/python_scripts/McSplicer.py \\
        --gtf ${refined_gtf} \\
        --count_file ${cnt_file} \\
        --out_dir output \\
        --read_len ${meta.read_len} \\
        ${gene_id_arg} \\
        ${gene_list_arg} \\
        ${prefix_arg} \\
        ${bootstrap_arg} \\
        ${anno_arg}
    """

    stub:
    """
    mkdir -p output/chr1
    touch output/chr1/test_gene.csv
    """
}
