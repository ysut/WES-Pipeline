nextflow.enable.dsl=2

params.output = "${projectDir}/../echtvar_zip"
params.prepro_data = "${projectDir}/processed_data"
params.chr_rename_rule = "${projectDir}/chr_rename_rules/chr2num.txt"
params.echtvar_json_37 = "${projectDir}/echtvar_json/gnomad4_remap37.json"
params.echtvar_json_38 = "${projectDir}/echtvar_json/gnomad4_GRCh38.json"


process EXTRACT_PASS_RENAME_CHR {
    publishDir params.prepro_data, mode: 'copy'

    input:
    tuple path(input_vcf_gz), path(index_file)

    output:
    path "*.pass.vcf.gz"

    script:
    base_name = file(input_vcf_gz).baseName.replace(".vcf", "")

    """
    bcftools view \\
      --threads 2 \\
      --output-type u \\
      --include 'FILTER="PASS"' \\
      ${input_vcf_gz} \\
      | bcftools annotate \\
          --threads 2 \\
          --rename-chrs ${params.chr_rename_rule} \\
          --output-type z \\
          > ${base_name}.pass.vcf.gz
    """
}

process CROSSMAP {
    input:
    tuple path(input_pass_vcf_gz), path(chain), path(target_ref)

    output:
    path '*.vcf'

    script:
    base_name = file(input_pass_vcf_gz).baseName.replace(".vcf", "")
    target_assembly = file(target_ref)
      .baseName
      .replace(".dna.primary_assembly.fa", "")
      .replace("Homo_sapiens.", "")

    """
    CrossMap vcf \
      ${chain} \
      ${input_pass_vcf_gz} \
      ${target_ref} \
      "${base_name}.${target_assembly}_remap.vcf"
    """
}

process EXCLUDE_NDot {
    publishDir params.prepro_data, mode: 'copy'
    
    input:
    path input_pass_remap_vcf

    output:
    path '*.nonNDot.bcf.gz'

    script:
    base_name = file(input_pass_remap_vcf).baseName.replace(".vcf", "")

    """
    bcftools view \\
      --threads 2 \\
      --output-type u \\
      --exclude 'REF="N" || REF="n" || ALT="N" || ALT="n" || ALT="."' \\
      ${input_pass_remap_vcf} \\
      | bcftools sort \\
          --max-mem 2G \\
          --output-type b \\
          > "${base_name}.nonNDot.bcf.gz"
    """
}

process INDEX_CSI {
    publishDir params.prepro_data, mode: 'move'

    input:
    path input_gz

    output:
    path '*.csi'

    script:
    """
    bcftools index --threads 2 --csi ${input_gz}
    """
}

process INDEX_CSI_FOR_REMAP {
    publishDir params.prepro_data, mode: 'move'

    input:
    path input_gz

    output:
    path '*.csi'

    script:
    """
    bcftools index --threads 2 --csi ${input_gz}
    """
}

process ENCODE_38 {
    publishDir params.output, mode: 'move'
    
    input:
    tuple path (input_data), path (config_json)

    output:
    path 'gnomad.v4.1_GRCh38.echtvar.zip'

    script:
    """
    echtvar encode \\
      gnomad.v4.1_GRCh38.echtvar.zip \\
      ${config_json} \\
      ${input_data}
    """
}

process ENCODE_37 {
    publishDir params.output, mode: 'move'
    
    input:
    tuple path (input_data), path (config_json)

    output:
    path 'gnomad.v4.1_GRCh37.echtvar.zip'

    script:
    """
    echtvar encode \\
      gnomad.v4.1_GRCh37.echtvar.zip \\
      ${config_json} \\
      ${input_data}
    """
}

workflow {
    // gnomad4 
    Channel.fromPath( "${params.input}/${params.gnomad_v4_1}" )
      | map { vcf_file ->
            def index_file = vcf_file + ".tbi"
            return [vcf_file, index_file]
        }
      | EXTRACT_PASS_RENAME_CHR
      | set { ch_pass_vcf_gz }

    // 1st. branch for non-remap data
    ch_pass_vcf_gz
      | INDEX_CSI

    ch_pass_vcf_gz
        .collect()
        .set { ch_collected_pass_vcf_gz }

    ch_collected_pass_vcf_gz
      | map { collected_pass_vcf_gz ->
                [collected_pass_vcf_gz, "${params.echtvar_json_38}"]
            }
      | ENCODE_38

    // 2nd. branch for remap data
    ch_pass_vcf_gz
      | map { pass_vcf_gz_pre_remap -> 
                [pass_vcf_gz_pre_remap, "${params.chain_to37}", "${params.ref_37}"]
            }
      | CROSSMAP
      | EXCLUDE_NDot
      | set { ch_remap_bcf_gz }

    ch_remap_bcf_gz  
      | INDEX_CSI_FOR_REMAP

    ch_remap_bcf_gz
        .collect()
        .set{ ch_collected_remap_bcf_gz }

    ch_collected_remap_bcf_gz
      | map { collected_pass_bcf_gz ->
                [collected_pass_bcf_gz, "${params.echtvar_json_37}"]
            }
      | ENCODE_37
}
