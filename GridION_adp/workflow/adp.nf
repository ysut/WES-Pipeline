#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_fastq = ''
params.trim_cutoff = 600

process TRIM_SHORT_READS {
    input:
    path input_fastq

    output:
    path 'trimmed.fastq'

    script:
    """
    seqkit seq -m ${params.trim_cutoff} ${input_fastq} > trimmed.fastq
    """
}

workflow {
    Channel.fromPath(params.input_fastq)
      | TRIM_SHORT_READS
      | view()
}