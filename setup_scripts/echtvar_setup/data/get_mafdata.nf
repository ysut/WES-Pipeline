params.output = ''
params.wes_dir = ''
params.fasta_release = 'release-112'

process DOWNLOADFASTA {
    publishDir "${params.output}/share", mode: 'move'

    base_uri = "https://ftp.ensembl.org/pub/${params.fasta_release}/fasta/homo_sapiens/dna"
    fasta_37 = "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
    fasta_38 = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    uri_37 = "${base_uri}/${fasta_37}"
    uri_38 = "${base_uri}/${fasta_38}"

    output:
    path 'Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz'
    path 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

    script:   
    """"
    curl -OkL ${uri_37} && gunzip ${fasta_37}
    curl -OkL ${uri_38} && gunzip ${fasta_38}
    bgzip -@ ${params.threads} Homo_sapiens.GRCh37.dna.primary_assembly.fa
    bgzip -@ ${params.threads} Homo_sapiens.GRCh38.dna.primary_assembly.fa
    """
}

process DOWNLOADCHAIN {
    publishDir "${params.output}/share", mode: 'move'

    chainfile_prefix = "https://hgdownload.soe.ucsc.edu/goldenPath"
    hg38ToHg19 = "${chainfile_prefix}/hg38/liftOver/hg38ToHg19.over.chain.gz"
    hg19ToHg38 = "${chainfile_prefix}/hg19/liftOver/hg19ToHg38.over.chain.gz"

    output:
    path 'hg38ToHg19.over.chain.gz'
    path 'hg19ToHg38.over.chain.gz'

    script:
    """
    


    """
}


process DOWNLOADCHRRULES {
    publishDir "${params.output}/share", mode: 'move'


    script:
    """
    """
}

process dbSNP {
    publishDir "${params.output}", mode: 'move'

    output:
    path 'GCF_000001405.25.dbSNP-156.singleALT.common.chrmod.sort.vcf.gz'
    path 'GCF_000001405.25.dbSNP-156.singleALT.common.chrmod.sort.vcf.gz.tbi'
    path 'GCF_000001405.40.dbSNP-156.singleALT.common.chrmod.sort.vcf.gz'
    path 'GCF_000001405.40.dbSNP-156.singleALT.common.chrmod.sort.vcf.gz.tbi'

    script:
    """
    ${params.wes_dir}/setup_workflow/ShellScripts/prep_dbSNP.sh
    """
}


process 


process RGCME {
    input:
    path fasta_37

    output:
    stdout

    script:
    """
    """
}

process EXCLUDEPASS {
    input:
    path 

    output:
    path 

    script:
    """
    bcftools view \\
      --threads 8 \\
      --output-type v \\
      --include 'FILTER="PASS"' \\
      > ${input_vcf_gz.basename}.vcf
    """
}





workflow {
    DOWNLOADFASTA()
        .set{fasta_37, fasta_38}

    RGCME(fasta_37)
    

}

