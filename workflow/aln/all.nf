#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.fastq_dir = ''
params.ped = ''
params.output = "${launchDir}/results"
// params.config_dir = "${projectDir}/config"
// params.vcfanno_conf_toml = "${params.config_dir}/vcfanno/conf_dbSNP_${params.ASSEMBLY}.toml"
// params.vep_options_ini = "${params.config_dir}/vep/options.ini"

// params.echtvar_gnomad4 = "${params.echtvar_resources}/gnomad.v4.1_${params.ASSEMBLY}.echtvar.zip"
// params.echtvar_tommo54k = "${params.echtvar_resources}/ToMMo54KJPN_auto_${params.ASSEMBLY}.echtvar.zip"
// params.echtvar_rgcme = "${params.echtvar_resources}/RGC_ME_${params.ASSEMBLY}.echtvar.zip"
// params.echtvar_ncbnwgs = "${params.echtvar_resources}/NCBN_${params.ASSEMBLY}.echtvar.zip"
// params.echtvar_genomeasia = "${params.echtvar_resources}/GenomeAsia_${params.ASSEMBLY}.echtvar.zip"

// params.phenotypes = ''

params.CHROMOSOME = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'
// params.CHROMOSOME = '18,22'


process PARSE_PED_FILE {
    publishDir "${params.output}/pedigree_files", mode: 'copy'

    input:
    path input_ped

    output:
    path '*.ped'

    script:
    """    
    while read -r line; do \\
      family_id=\$(echo "\$line" | awk '{print \$1}') \\
      sanitized_family_id=\$(echo "\${family_id}" | sed 's/[\\/|:*?"<> ]/_/g') \\
      echo "\$line" >> "\${sanitized_family_id}.ped" \\
    done < "${input_ped}"
    """
}

process FIND_FASTQ_FILES {
    input:
    path '*.ped'

    output:
    tuple val(sample_id), path(read1), path(read2)

    script:
    """
    while read sample_id; do
        read1=\$(find ${params.fastq_dir} -type f -name "*\${sample_id}*_R1*.fastq.gz")
        read2=\$(find ${params.fastq_dir} -type f -name "*\${sample_id}*_R2*.fastq.gz")
        echo -e "${sample_id}\\t${read1}\\t${read2}" >> found_fastqs.txt
    done < *.ped
    """
}


// process STROBEALIGN {
//     // container 'betelgeuse:5000/library/utsu/strobealign:0.14.0'
    
//     input:
//     tuple path(read1), path(read2), path(reference), val(sample_id)
    
//     output:
//     tuple path('sorted.bam'), path('sorted.bam.bai'), val(sample_id)

//     script:
//     """
//     conda activate strobealign && \\
    
//     strobealign \\
//       --threads=8 ${reference} ${read1} ${read2} 
//       --rg-id=${sample_id} --rg=SM:${sample_id} --rg=LB:mylibrary --rg=PL:Illumina
//         | samtools sort -o sorted.bam && \\
    
//     samtools index sorted.bam
//     """
// }

process MARKDUP {
    publishDir "${params.output}/markdup_files", mode: 'copy', pattern: '*.txt'
    // container 'betelgeuse:5000/library/utsu/samtools:1.10'
    // docker run --rm e0759ef897a1 java -jar /usr/picard/picard.jar MarkDuplicates --help

    input:
    tuple path(aln_file), path(aln_index), val(sample_id)
    
    output:
    tuple path('sorted.markdup.bam'), path('sorted.markdup.bam.bai'), val(sample_id)

    script:
    """
    java -jar /usr/picard/picard.jar MarkDuplicates \\
      -I ${aln_file} \\
      -O sorted.markdup.bam \\
      -M ${sample_id}_marked_dup_metrics.txt \\
    samtools index sorted.markdup.bam
    """
}


process DEEPVARIANT {
    container 'betelgeuse:5000/library/utsu/deepvariant:1.6.0'
    containerOptions '-v /betelgeuse07/analysis/utsu/resources:/resources'

    input:
    tuple(val(sample_id), path(aln_file), path(aln_index))

    output:
    tuple(path("${sample_id}.g.vcf.gz"), path("${sample_id}.g.vcf.gz.tbi"))

    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref=${fasta} \
    --reads=${aln_file} \
    --output_vcf="${sample_id}.vcf.gz" \
    --output_gvcf="${sample_id}.g.vcf.gz" \
    --num_shards=64 \
    --intermediate_results_dir=intermediate \
    --logging_dir=logs
    """
}



workflow {
    Channel.fromPath(params.ped)
      | PARSE_PED_FILE
      | FIND_FASTQ_FILES
      | view()
      // | map { it -> tuple(it, "${params.fasta}") }
    // // Using split files by chromosome
    // ch_vcfs
    //   | map {it -> tuple(it, "${params.fasta}", "${params.spliceai_annotation}") }
    //   | SPLICEAI
    //   | set {ch_spliceai_vcfs}

    // ch_spliceai_vcfs.collect()
    //   | CONCATENATE_VCFS
    //   | map { vcfs -> tuple(vcfs, "${params.fasta}", "${params.vep_options_ini}") }
    //   | VEP
    //   | set{ vep_output_ch }
}

