#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input = ''
params.ped = ''
params.output = "${launchDir}/results"
params.config_dir = "${projectDir}/config"
params.vcfanno_conf_toml = "${params.config_dir}/vcfanno/conf_dbSNP_${params.ASSEMBLY}.toml"
params.vep_options_ini = "${params.config_dir}/vep/options.ini"

params.echtvar_gnomad4 = "${params.echtvar_resources}/gnomad.v4.1_${params.ASSEMBLY}.echtvar.zip"
params.echtvar_tommo54k = "${params.echtvar_resources}/ToMMo54KJPN_auto_${params.ASSEMBLY}.echtvar.zip"
params.echtvar_rgcme = "${params.echtvar_resources}/RGC_ME_${params.ASSEMBLY}.echtvar.zip"
params.echtvar_ncbnwgs = "${params.echtvar_resources}/NCBN_${params.ASSEMBLY}.echtvar.zip"
params.echtvar_genomeasia = "${params.echtvar_resources}/GenomeAsia_${params.ASSEMBLY}.echtvar.zip"

// params.phenotypes = ''

params.CHROMOSOME = '1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,MT'
// params.CHROMOSOME = '18,22'

process PARSE_PED_FILE {
    input:
    path input_ped

    output:
    // tupleかリストで，famliy_id，proband_id，father_id，mother_id返す
    // Ifで制御必要，

}


process EXTRACT_CHR {
    // Extract chromosome from the input VCF file
    input:
    tuple path(input_vcf), path(fasta)

    output:
    path 'chr.vcf'
      // --fasta-ref /betelgeuse05/system/bio/db/riker/bundle_b37/human_g1k_v37_fix.fasta \\

    script:
    // In nextflow, bcftools can accept only one input file.
    """
    bgzip \\
      --threads ${params.threads} \\
      ${input_vcf}

    tabix \\
      --preset vcf \\
      ${input_vcf}.gz
    
    bcftools view \\
      --threads ${params.threads} \\
      --regions ${params.CHROMOSOME} \\
      --output-type b \\
      ${input_vcf}.gz \\
      > "chr.vcf"      
    """
}

process NORM {
    input:
    path input_bcf

    output:
    path 'chr.norm.vcf.gz'

    script:
    """
    bcftools norm \\
      --multiallelics -both \\
      --threads ${params.threads} \\
      --fasta-ref ${params.fasta} \\
      --output-type z \\
      ${input_bcf} \\
      > chr.norm.vcf.gz
    """
}

process VCFANNO_DBSNP {
    input:
    path input_vcf

    output:
    path 'chr.dbsnp.vcf'

    script:
    """
    vcfanno \\
      -p ${params.vcfanno_threads} \\
      ${params.vcfanno_conf_toml} \\
      ${input_vcf} \\
      > chr.dbsnp.vcf
    """
}

process EXCLUDE_COMMON {
    input:
    path input_vcf

    output:
    path 'chr.uncommon.bcf'

    script:
    """
    bcftools view \\
      --threads ${params.threads} \\
      --exclude 'dbSNPcommon=1' \\
      --output-type u \\
      ${input_vcf} \\
      > chr.uncommon.bcf \\
    """
}

process ECHTVAR_ANNO {
    // Annotate MAF data with echtvar
    input:
    path input_bcf

    output:
    path 'chr.uncommon.echtvar.bcf'

    script:
    """
    echtvar anno \\
      -e ${params.echtvar_gnomad4} \\
      -e ${params.echtvar_tommo54k} \\
      -e ${params.echtvar_rgcme} \\
      -e ${params.echtvar_ncbnwgs} \\
      -e ${params.echtvar_genomeasia} \\
      ${input_bcf} \\
      chr.uncommon.echtvar.bcf
    """
}


process VCFANNO_HGMD {    
    input:
    path input_vcf

    output:
    path 'chr.uncommon.echtvar.hgmd.bcf'

    script:
    """
    vcfanno \\
      -p ${params.vcfanno_threads} \\
      ${params.vcfanno_conf_toml} \\
      ${input_vcf} \\
      > chr.uncommon.hgmd.bcf
    """
}

process MAF_FILTER {
    // Extract rare variants using bcftools by MAF cutoff
    input:
    path input_bcf

    output:
    tuple path('chr.uncommon.rare.vcf.gz'), path('chr.uncommon.rare.vcf.gz.csi')

    script:
    """
    bcftools view \\
      --threads ${params.threads} \\
      --include 'INFO/gnomad4_remap37_af_grpmax_joint<=${params.MAFCUTOFF} && INFO/tommo54kjpn_lifted37_af<=${params.MAFCUTOFF} && INFO/rgcme_remap37_af_all<=${params.MAFCUTOFF} && INFO/ncbn_remap37_af_jpn<=${params.MAFCUTOFF} && INFO/ga_af<=${params.MAFCUTOFF}' \\
      --output-type z \\
      ${input_bcf} \\
      > chr.uncommon.rare.vcf.gz

    bcftools index \\
      --csi \\
      --threads ${params.threads} \\
      chr.uncommon.rare.vcf.gz
    """
}

process SPLIT_BY_CHR {
    input:
    path vcf

    output:
    path 'chr.*.vcf'

    script:
    """
    for chr in {1..22} X Y MT; do
      bcftools view \\
        --threads ${params.threads} \\
        --regions \${chr} \\
        --output-type v \\
        ${vcf} \\
        > chr.\${chr}.vcf
    done
    """
}

process SPLICEAI {
    input:
    tuple path(input_vcf), path(fasta), path(annotation)

    output:
    path '*.splai.vcf'

    script:
    """
    source /opt/conda/etc/profile.d/conda.sh && \\
    conda activate spliceai && \\
    spliceai \\
      -I ${input_vcf} \\
      -O ${input_vcf.baseName}.splai.vcf \\
      -R ${fasta} \\
      -A ${annotation} \\
      -D ${params.spliceai_distance} \\
      -M ${params.spliceai_mask}
    """
}

process CONCATENATE_VCFS {
    input:
    path collected_vcfs

    output:
    path 'chr.uncommon.rare.spliceai.vcf'

    script:
    """
    which bcftools
    bcftools concat \\
      --threads ${params.threads} \\
      --output-type v \\
      ${collected_vcfs} \\
      > chr.uncommon.rare.spliceai.vcf \\
    """
}

process VEP {
    // publishDir "/home/utsu/Desktop", mode: 'copy'
    publishDir "${params.output}", mode: 'copy'
    containerOptions "-u 0 -v ${params.vep_data}:/data -v ${params.vep_plugin_resources}:/plugin_resources"

    input:
    tuple path(input_vcf), path(fasta), path(options_ini)

    output:
    // tuple(path("nf.vep.txt"), path("nf.vep.txt_summary.html"))
    path 'nf.vep.vcf'

    script:
    """ 
    /opt/vep/src/ensembl-vep/vep \\
      --input_file ${input_vcf} \\
      --output_file nf.vep.vcf \\
      --fasta ${fasta} \\
      --config ${options_ini} \\
      --assembly ${params.ASSEMBLY} \\
      --buffer_size ${params.vep_buffer} \\
      --fork ${params.threads}
    """
}

process PARSE_XHMM {
    input:
    path xhmm_output
    // var proband_id を受け取って，joint callされているファイルからproband_idが
    // 一致する行だけ抽出してproband_id.xhmm.txt　みたいのを作る
    // outputで返す．あとでmergeする．

    output:
    path 'xhmm_summary.txt'

    script:
    """
    grep -A 1 'Z-score' ${xhmm_output} \\
      | tail -n 1 \\
      > xhmm_summary.txt
    """
}



process JIGV {
    publishDir "${params.output}/", mode: 'move'
    input:
    tuple path(input_vcf), path(fasta), path(ped), path(vep_vcf)

    script:
    """
    jigv \\
      --sample ${sample} \\
      --ped ${ped} \\
      --sites 
      --fasta ${fasta} \\
    echo "hoge"
    """
}

workflow.onComplete {
      println ""
      println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
      println "Pipeline completed at: $workflow.complete"
      println "Execution time       : $workflow.duration"
      println "Execution status     : ${ workflow.success ? 'OK' : 'failed' }"
      println "~*~*~*~*~~*~*~*~*~~*~*~*~*~~*~*~*~*~*~*~*~~*~*~*~*~*~*~*~*~"
}

workflow {
    Channel.fromPath(params.input)
      | map { input_vcf -> [input_vcf, "${params.fasta}"] }
      // | map { input_vcf -> tuple(input_vcf, "${params.fasta}", "${params.fasta_index}") }
      | EXTRACT_CHR
      | NORM
      | VCFANNO_DBSNP
      | EXCLUDE_COMMON
      | ECHTVAR_ANNO
      | MAF_FILTER
      | SPLIT_BY_CHR
      | flatMap { vcf -> vcf }
      | set {ch_vcfs}
    
    // Using split files by chromosome
    ch_vcfs
      | map {it -> tuple(it, "${params.fasta}", "${params.spliceai_annotation}") }
      | SPLICEAI
      | set {ch_spliceai_vcfs}

    ch_spliceai_vcfs.collect()
      | CONCATENATE_VCFS
      | map { vcfs -> tuple(vcfs, "${params.fasta}", "${params.vep_options_ini}") }
      | VEP
      | set{ vep_output_ch }
   
}
