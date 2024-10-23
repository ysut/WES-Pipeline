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

process FIND_FASTQ_FILES {
  input:
    val sample_id 

  output:
    tuple val(sample_id), path(read1), path(read2)

  script:
  """
  # サンプルIDを使用して対応するFASTQファイルを探す
  read1=\$(find . -name "${sample_id}*R1*.fastq.gz")
  read2=\$(find . -name "${sample_id}*R2*.fastq.gz")

  # 見つかったパスを返す
  echo "$sample_id \$read1 \$read2" > fastq_files.txt
  """
}


workflow {
    // Load a pedigree format file
    ch_sample_ids = Channel.fromPath(params.ped)
      .splitText()
      .map { line -> line.split('\t')[1] } // 2列目のサンプルIDを抽出
      // .set { sample_ids }
    
    FIND_FASTQ_FILES(ch_sample_ids)
      | view()
    
}

