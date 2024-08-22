#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// params.input = ''
// params.phenotypes = ''

// The path below is container specific 
params.fasta = '/resources/ReferenceGenome/exome_pipeline/human_g1k_v37_fix.fasta'

process GETSAMPLES {
    // Find all bam/bai files in the "preparation" directory
    // ToDo: cram and crai (bam.bai/caram.crai) files should be included
    output:
    path "samples.txt"

    script:
    """
    work_dir=${workflow.workDir}
    aln_root_dir=\${work_dir%/*}/preparation
    find \$aln_root_dir -type f -name "*.bam" > alns.txt
    
    while read aln_path; do
      basename=\$(basename \${aln_path})
      sample_id=\${basename%%.*}
      index_path=\$(echo \${aln_path} | sed 's/am\$/ai/')
      echo "\${sample_id}\t\${aln_path}\t\${index_path}" >> samples.txt
    done < alns.txt
    """
}

process DEEPVARIANT {
    container 'betelgeuse:5000/library/utsu/deepvariant:1.6.0'
    containerOptions '-v /betelgeuse07/analysis/utsu/resources:/resources'

    input:
        tuple(val(sample_id), path(aln_file), path(aln_index))
        val fasta

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

process GLNEXUS {
    container 'betelgeuse:5000/library/utsu/glnexus:1.4.3'

    input:
        path(gvcf_files)
        path(gvcf_indices)
    
    output:
        path "joint.bcf"
    
    script:
    """
    /bin/bash -c 'glnexus_cli \
    --config DeepVariantWES \
    --mem-gbytes 256 \
    --threads 64 \
    ${gvcf_files}' > joint.bcf
    """
}

process BCF2VCF {
    input:
        path bcf_file

    output:
        tuple(path("joint.vcf.gz"), path("joint.vcf.gz.tbi"))

    script:
    """
    bcftools view \\
      --output-type z \\
      --output-file joint.vcf.gz \\
      --threads 16 \\
      ${bcf_file}

    bcftools index \\
      --tbi \\
      --threads 16 \\
      joint.vcf.gz
    """
}

process VCFANNO {
    input:
    tuple(path(vcf_file), path(vcf_index))

    output:
    path 'vcfanno.vcf'

    script:
    """
    vcfanno \\
      -p 8 \\
      ${WORKFLOW_WES}/config/vcfanno_${ASSEMBLY}.toml \\
      $vcf_file \\
      > vcfanno.vcf
    """
}

process SPLIT {
    input:
    path vcf

    output:
    path 'chunk_*'

    script:
    """
    ${WORKFLOW_WES}/scripts/chrsplit.sh $vcf
    """
}

process JPN575 {
    input:
        path vcf
    
    output:
        path "*.inhouse.vcf"

    script:
    """   
    perl /usr/local/bits/riker3/bin/japanese_mutation_filter.pl \\
      $vcf > ${vcf.baseName}.inhouse.vcf
    """
}

process SPLICEAI {
    container 'betelgeuse:5000/library/utsu:spliceai'
    containerOptions '-v /betelgeuse07/analysis/utsu/resources:/resources'
    
    input:
        path vcf

    output:
        path '*.splai.vcf'

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh
    conda activate spliceai
    spliceai \\
      -I $vcf \\
      -O ${vcf.baseName}.splai.vcf \\
      -R /resources/ReferenceGenome/exome_pipeline/human_g1k_v37_fix.fasta \\
      -A /resources/SAI10kcalc/SAI10kcalc_spliceai_tx.hg19.tsv \\
      -D 4999 \\
      -M 0
    """
}

process CONCATENATESORT {
    input:
    path collected_vcfs

    output:
    path 'splai.concat.sort.pre_annovar.vcf'

    script:
    """
    bcftools concat \\
      --output splai.concat.pre_annovar.vcf \\
      --output-type v \\
      --threads 4 \\
      ${collected_vcfs}

    bcftools sort \\
      --output-file splai.concat.sort.pre_annovar.vcf \\
      --output-type v \\
      --max-mem 16G \\
      splai.concat.pre_annovar.vcf
    """
}


process ANNOVAR {
    input:
    path pre_annovar_vcf

    output:
    tuple path( '*.txt' ), path( '*.vcf' ), emit: annovar_files

    script:
    """
    ${WORKFLOW_WES}/scripts/run_annovar_wes_${ASSEMBLY}.sh \\
      ${pre_annovar_vcf} \\
      ${ANNOVAR_DB} \\
      ${SPLICING_THRESHOLD}
    """
}

process FORMATANNOVAR {
    input:
    tuple path(txt), path(vcf)

    output:
    path "${txt.baseName}.renamed.txt"

    script:
    """
    ${WORKFLOW_WES}/bin/anvrformatter \\
      ${txt} \\
      ${WORKFLOW_WES}/config/rename.toml
    """
}

process HGMDANNOTATOR {
    input:
    path txt

    output:
    path 'exome_summary_*.txt'

    script:
    """
    run_hgmd_annotator ${txt} exome_summary.txt
    DATE=\$(date +'%Y%m%d_%H%M%S')
    mv exome_summary.txt exome_summary_\${DATE}.txt
    cp exome_summary_\${DATE}.txt ${workflow.launchDir}/
    """
}

process GADO {
    input:
    path phenotypes

    output:
    path 'GADO_results/*.txt'

    script:
    """
    ${WORKFLOW_WES}/scripts/run_gado.sh $phenotypes ./GADO_results/
    mkdir ./GADO_results/log
    mv ./GADO_results/hpoProcessed.txt* ./GADO_results/log/
    mv ./GADO_results/samples.txt ./GADO_results/log/
    cp ./GADO_results/*.txt ${workflow.launchDir}/
    """
}

workflow {
    fasta = "${file(params.fasta).toAbsolutePath()}"

    GETSAMPLES()
        .set{ samples_ch }

    samples_ch
        .splitCsv(sep: '\t', strip: true)
        .map { it -> tuple(it[0], it[1], it[2])}
        .set { sample_info_ch }
    
    DEEPVARIANT(sample_info_ch, fasta)
        .set{ deepvariantFiles_ch }

    deepvariantFiles_ch
        .map{ it[0] }
        .collect()
        .set{ gvcf_files_ch }

    deepvariantFiles_ch
        .map{ it[1] }
        .collect()
        .set{ gvcf_indices_ch }

    GLNEXUS(gvcf_files_ch, gvcf_indices_ch)
        .set{ glnexusFile_ch }

    BCF2VCF(glnexusFile_ch)
        .set{ vcfFile_ch }
    
    VCFANNO(vcfFile_ch)
        .set{ vcfannoFile_ch }

    SPLIT(vcfannoFile_ch)
        .flatMap{ file -> file }
        .set{ splitFiles_ch }

    JPN575(splitFiles_ch)
        .set{ jpn575File_ch }

    SPLICEAI(jpn575File_ch)
        .collect()
        .set{ spliceaiFiles_ch }

    spliceaiFiles_ch
        .collect()
        .set{ collectedSpliceaiFiles_ch }
    
    CONCATENATESORT(collectedSpliceaiFiles_ch)
        .set{ concatFile_ch }

    if (params.phenotypes) {
        phenotypes_ch = Channel.fromPath(params.phenotypes)
        GADO(phenotypes_ch)
            .set{ gadoFile_ch }
    }
    
    // MAVERICK(input_ch)
    //     .set{ maverickFiles_ch }
    // maverickFiles_ch.view()

    ANNOVAR(concatFile_ch)
        .set{ annovarFiles_ch }

    FORMATANNOVAR(annovarFiles_ch)
        .set{ renamedFile_ch }

    HGMDANNOTATOR(renamedFile_ch)
        .set{ hgmdAnnotatedFile_ch }

    hgmdAnnotatedFile_ch.view()    

}
