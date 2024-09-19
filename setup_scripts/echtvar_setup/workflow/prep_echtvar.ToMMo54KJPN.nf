// ToMMo 54KJPN were build on GRCh38 and GRCh37 lifted from GRCh38
// Composed data. Without chr prefix.

nextflow.enable.dsl=2

params.output = "${projectDir}/../echtvar_zip"
params.prepro_data = "${projectDir}/processed_data"
params.echtvar_json_37 = "${projectDir}/echtvar_json/tommo54kjpn_auto_lifted37.json"
params.echtvar_json_38 = "${projectDir}/echtvar_json/tommo54kjpn_auto_GRCh38.json"


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

process EXTRACT_PASS_NORM_SORT_FOR_AUTO37 {
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
          --multiallelics -both \\
          --output-type u \\
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

process INDEX_CSI_FOR_AUTO38 {
    publishDir params.prepro_data, mode: 'move'

    input:
    path input_gz

    output:
    path '*.csi'

    script:
    """
    bcftools index --csi ${input_gz}
    """
}

process INDEX_CSI_FOR_AUTO37 {
    publishDir params.prepro_data, mode: 'move'

    input:
    path input_gz

    output:
    path '*.csi'

    script:
    """
    bcftools index --csi ${input_gz}
    """
}

process ENCODE_37 {
    publishDir params.output, mode: 'move'
    
    input:
    tuple path (input_data), path (config_json)

    output:
    path 'ToMMo54KJPN_auto_GRCh37.echtvar.zip'

    script:
    """
    echtvar encode \\
      ToMMo54KJPN_auto_GRCh37.echtvar.zip \\
      ${config_json} \\
      ${input_data}
    """
}

process ENCODE_38 {
    publishDir params.output, mode: 'move'
    
    input:
    tuple path (input_data), path (config_json)

    output:
    path 'ToMMo54KJPN_auto_GRCh38.echtvar.zip'

    script:
    """
    echtvar encode \\
      ToMMo54KJPN_auto_GRCh38.echtvar.zip \\
      ${config_json} \\
      ${input_data}
    """
}


workflow {
    // ToMMo 54KJPN autosomal on GRCh38
    Channel.fromPath( "${params.input}/${params.tommo54kjpn_38_auto}" )
      | map { vcf_file ->
            def index_file = vcf_file + ".tbi"
            return [vcf_file, index_file]
        }
      | EXTRACT_PASS_NORM_SORT_FOR_AUTO38
      | set {ch_pass_38_bcf_gz}

    ch_pass_38_bcf_gz
      | INDEX_CSI_FOR_AUTO38
    
    ch_pass_38_bcf_gz
      | map { pass_bcf_gz ->
                [pass_bcf_gz, "${params.echtvar_json_38}"]
            }
      | ENCODE_38

    // ToMMo 54KJPN autosomal on GRCh37
    Channel.fromPath( "${params.input}/${params.tommo54kjpn_37_auto}" )
      | map { vcf_file ->
            def index_file = vcf_file + ".tbi"
            return [vcf_file, index_file]
        }
      | EXTRACT_PASS_NORM_SORT_FOR_AUTO37
      | set {ch_pass_37_bcf_gz}

    ch_pass_37_bcf_gz
      | INDEX_CSI_FOR_AUTO37
    
    ch_pass_37_bcf_gz
      | map { pass_bcf_gz ->
                [pass_bcf_gz, "${params.echtvar_json_37}"]
            }
      | ENCODE_37
}
