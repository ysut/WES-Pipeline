#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.input_fastq = ''
params.input_trimmed_fastq = ''
params.trim_cutoff = 550

params.reference_fasta = '/betelgeuse07/analysis/utsu/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa'

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
    tuple path('aligned.sorted.bam'), path('aligned.sorted.bam.bai')

    script:
    """
    minimap2 -ax map-ont -t 12 ${reference_fasta} ${input_fastq} \\
      | samtools view -b - \\
      | samtools sort -m 4G -@ 8 -O bam -o aligned.sorted.bam - && \\
    samtools index aligned.sorted.bam
    """
}

process LONGSHOT {
    publishDir "${launchDir}", mode: 'move'

    input:
    tuple path(input_bam), path(input_bai). path(reference_fasta), path(reference_fasta_index)

    output:
    path 'aligned.sorted.longshot.bam'

    script:
    """
    source /usr/local/genome/miniconda3/miniconda/bin/activate longshot-0.4.1 && \\
    longshot \\
      --auto_max_cov \\
      --bam ${input_bam} \\
      --ref ${reference_fasta} \\
      --strand_bias_pvalue_cutoff 0.01 \\
      --out longshot.vcf \\
      --out_bam aligned.sorted.longshot.bam
    """
}

workflow {
    // Channel.fromPath(params.input_fastq)
    //   | TRIM_SHORT_READS
    Channel.fromPath(params.input_trimmed_fastq)
      | map {it -> tuple(it, "${params.reference_fasta}")}
      | MINIMAP2
      | map {it -> tuple(it[0], it[1], "${params.reference_fasta}", "${params.reference_fasta}.fai")}
      | LONGSHOT
      | view()
}