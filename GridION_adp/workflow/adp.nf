#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_fastq = ''
params.input_trimmed_fastq = ''
params.trim_cutoff = 500

params.reference_fasta = '/betelgeuse01/analysis/miyatake/minimap2/GRCh38_full_analysis_set_plus_decoy_hla.fa'

process TRIM_SHORT_READS {
    publishDir "${launchDir}", mode: 'move'

    input:
    path input_fastq

    output:
    path 'trimmed.fastq'

    script:
    """
    seqkit seq -m ${params.trim_cutoff} --threads 8 ${input_fastq} > trimmed.fastq
    """
}

process MINIMAP2 {
    input:
    tuple path(input_fastq), path(reference_fasta)

    output:
    path 'aligned.bam'

    script:
    """
    minimap2 -ax map-ont -t 12 ${reference_fasta} ${input_fastq} \\
      | samtools view -b - \\
      | samtools sort -m 1G -@ 8 -O bam -o aligned.bam -
    """
}

workflow {
    // Channel.fromPath(params.input_fastq)
    //   | TRIM_SHORT_READS
    Channel.fromPath(params.input_trimmed_fastq)
      | map {it -> tuple(it, "${params.reference_fasta}")}
      | MINIMAP2
      | view()
}