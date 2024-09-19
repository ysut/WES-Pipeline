// ToMMo 54KJPN were build on GRCh38 and GRCh37 lifted from GRCh38
// Composed data. Without chr prefix.

nextflow.enable.dsl=2

params.output = "${projectDir}/../echtvar_zip"
params.prepro_data = "${projectDir}/processed_data"
params.echtvar_json_37 = "${projectDir}/echtvar_json/hgmd_hg19.json"
params.echtvar_json_38 = "${projectDir}/echtvar_json/hgmd_hg38.json"


process EXTRACT_PASS_NORM_SORT_FOR_AUTO38 {
    publishDir params.prepro_data, mode: 'copy'

    input:
    tuple path(input_vcf_gz), path(index_file)

    output:
    path "*.pass.singleALT.nonNDot.bcf.gz"

    script:
    base_name = file(input_vcf_gz).baseName.replace(".vcf", "")

    """
    bcftools view \\
      --threads 4 \\
      --output-type u \\
      --include 'FILTER="PASS"' \\
      ${input_vcf_gz} \\
      | bcftools norm \\
          --threads 4 \\
          --output-type u \\
          --multiallelics -both \\
      | bcftools view \\
          --threads 4 \\
          --output-type u \\
          --exclude 'REF="N" || REF="n" || ALT="N" || ALT="n" || ALT="."' \\
      | bcftools sort \\
          --max-mem ${params.sort_max_mem} \\
          --output-type b \\
          > "${base_name}.pass.singleALT.nonNDot.bcf.gz"
    """
}


process ENCODE_37 {
    publishDir params.output, mode: 'move'
    
    input:
    tuple path (input_data), path (config_json)

    output:
    path 'HGMD_*_hg19.echtvar.zip'

    script:
    """
    echtvar encode \\
      HGMD_${version}_hg19.echtvar.zip \\
      ${config_json} \\
      ${input_data}
    """
}

process ENCODE_38 {
    publishDir params.output, mode: 'move'
    
    input:
    tuple path(input_data), path(config_json)

    output:
    path 'HGMD_*_hg38.echtvar.zip'

    script:
    """
    echtvar encode \\
      HGMD_${version}_hg38.echtvar.zip \\
      ${config_json} \\
      ${input_data}
    """
}


workflow {
    // hg38
    Channel.fromPath( "${params.input}/${params.tommo54kjpn_38_auto}" )
      | map { pass_vcf_gz ->
                [pass_vcf_gz, "${params.echtvar_json_38}"]
            }
      | ENCODE_38

    // hg19
    Channel.fromPath( "${params.input}/${params.tommo54kjpn_37_auto}" )
      | map { pass_vcf_gz ->
                [pass_vcf_gz, "${params.echtvar_json_37}"]
            }
      | ENCODE_37
}
