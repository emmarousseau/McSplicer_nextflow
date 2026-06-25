#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    McSplicer — Splice Site Usage Estimation Pipeline
========================================================================================
    Steps: exonRefine -> sigcount -> McSplicer EM
    https://github.com/mcsplicer
========================================================================================
*/

include { EXON_REFINE } from './modules/local/exon_refine'
include { SIGCOUNT    } from './modules/local/sigcount'
include { MCSPLICER   } from './modules/local/mcsplicer'

/*
========================================================================================
    PARAMETERS
========================================================================================
*/

params.input            = null
params.outdir           = './results'
params.read_len         = 151
params.gene_id          = ''
params.prefix           = ''
params.bootstraps       = 0
params.anno             = 'n'
params.paired_end       = false
params.skip_preprocess  = false

/*
========================================================================================
    WORKFLOW
========================================================================================
*/

workflow {

    if (!params.input) {
        error "ERROR: Please provide a samplesheet via --input <path/to/samplesheet.tsv>"
    }

    // Parse samplesheet
    ch_input = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def meta = [:]
            meta.id              = row.sample_id
            meta.read_len        = row.read_len        ? row.read_len.toInteger()    : params.read_len
            meta.gene_id         = row.gene_id         ?: params.gene_id
            meta.prefix          = row.prefix          ?: params.prefix
            meta.paired_end      = row.paired_end      ? row.paired_end.toLowerCase() in ['true','yes','1'] : params.paired_end
            meta.skip_preprocess = row.skip_preprocess ? row.skip_preprocess.toLowerCase() in ['true','yes','1'] : params.skip_preprocess
            meta.bootstraps      = row.bootstraps      ? row.bootstraps.toInteger()  : params.bootstraps
            meta.anno            = row.anno            ?: params.anno

            def bam = file(row.bam_file, checkIfExists: true)
            def bai = file("${row.bam_file}.bai", checkIfExists: true)
            def gtf = file(row.gtf_file, checkIfExists: true)

            return [ meta, bam, bai, gtf ]
        }

    // Branch: samples that need preprocessing vs those that skip it
    ch_input.branch {
        preprocess: !it[0].skip_preprocess
        skip:        it[0].skip_preprocess
    }.set { ch_branched }

    // ──────────────────────────────────────────────────────────────────────────
    // STEP 1: exonRefine — deduplicated per unique GTF
    // ──────────────────────────────────────────────────────────────────────────

    ch_unique_gtfs = ch_branched.preprocess
        .map { meta, bam, bai, gtf -> gtf }
        .unique()

    EXON_REFINE( ch_unique_gtfs )

    // ──────────────────────────────────────────────────────────────────────────
    // STEP 2: sigcount — per sample
    // ──────────────────────────────────────────────────────────────────────────

    // Join refined GTF back to each sample by matching original GTF path
    ch_for_sigcount = ch_branched.preprocess
        .map { meta, bam, bai, gtf -> [ gtf.name, meta, bam, bai ] }
        .combine(
            EXON_REFINE.out.refined_gtf.map { orig, refined -> [ orig.name, refined ] },
            by: 0
        )
        .map { gtf_name, meta, bam, bai, refined_gtf -> [ meta, bam, bai, refined_gtf ] }

    SIGCOUNT( ch_for_sigcount )

    // ──────────────────────────────────────────────────────────────────────────
    // STEP 3: McSplicer EM — per sample
    // ──────────────────────────────────────────────────────────────────────────

    // Preprocessed samples
    ch_mcsplicer_preproc = SIGCOUNT.out.counts
        .map { meta, cnt_file, refined_gtf -> [ meta, refined_gtf, cnt_file ] }

    // Skipped samples: GTF is already refined, .cnt is co-located with BAM
    ch_mcsplicer_skip = ch_branched.skip
        .map { meta, bam, bai, gtf ->
            def cnt_file = file(bam.toString().replaceAll(/\.bam$/, '.cnt'), checkIfExists: true)
            [ meta, gtf, cnt_file ]
        }

    ch_for_mcsplicer = ch_mcsplicer_preproc.mix(ch_mcsplicer_skip)

    MCSPLICER( ch_for_mcsplicer )
}
