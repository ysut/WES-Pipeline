params.output = ''
params.primateai_dir = ''
params.dbnsfp_version = '4.9a'

process LOFTEE {
    publishDir "${params.output}", mode: 'move'
    
    output:
    path 'human_ancestor.fa.gz'
    path 'human_ancestor.fa.gz.fai'
    path 'human_ancestor.fa.gz.gzi'
    path 'phylocsf_gerp.sql.gz'

    script:
    """
    curl -OkL https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz
    curl -OkL https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.fai
    curl -OkL https://s3.amazonaws.com/bcbio_nextgen/human_ancestor.fa.gz.gzi
    curl -OkL https://personal.broadinstitute.org/konradk/loftee_data/GRCh37/phylocsf_gerp.sql.gz

    """
}



process REVEL {
    publishDir "${params.output}", mode: 'move'
    
    output:
    path 'new_tabbed_revel_GRCh37.tsv.gz'
    path 'new_tabbed_revel_GRCh37.tsv.gz.tbi'
    path 'new_tabbed_revel_GRCh38.tsv.gz'
    path 'new_tabbed_revel_GRCh38.tsv.gz.tbi'

    script:
    """
    curl -OkL# https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip
    unzip revel-v1.3_all_chromosomes.zip
    cat revel_with_transcript_ids \\
      | tr "," "\t" > tabbed_revel.tsv
    sed '1s/.*/#&/' tabbed_revel.tsv > new_tabbed_revel_GRCh37.tsv
    bgzip new_tabbed_revel_GRCh37.tsv
    tabix -f -s 1 -b 2 -e 2 new_tabbed_revel_GRCh37.tsv.gz
    zcat new_tabbed_revel_GRCh37.tsv.gz | head -n1 > h
    zgrep -h -v ^#chr new_tabbed_revel_GRCh37.tsv.gz \\
      | awk '\$3 != "." ' \\
      | sort -k1,1 -k3,3n - \\
      | cat h - \\
      | bgzip -c > new_tabbed_revel_GRCh38.tsv.gz
    tabix -f -s 1 -b 3 -e 3 new_tabbed_revel_GRCh38.tsv.gz
    """
}

process ALPHAMISSENSE {
    publishDir "${params.output}", mode: 'move'

    output:
    path 'AlphaMissense_hg19.tsv.gz'
    path 'AlphaMissense_hg19.tsv.gz.tbi'
    path 'AlphaMissense_hg38.tsv.gz'
    path 'AlphaMissense_hg38.tsv.gz.tbi'

    script:
    """
    curl -OkL https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg19.tsv.gz
    curl -OkL https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
    tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg38.tsv.gz
    tabix -s 1 -b 2 -e 2 -f -S 1 AlphaMissense_hg19.tsv.gz
    """
}

process DOSAGESENSITIVITY {
    publishDir "${params.output}", mode: 'move'

    output:
    path 'Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz'

    script:
    """
    curl -OkL https://zenodo.org/record/6347673/files/Collins_rCNV_2022.dosage_sensitivity_scores.tsv.gz
    """
}

process EMFORMER {
    script:
    """
    curl -OkL https://ftp.ensembl.org/pub/current_variation/Enformer/enformer_grch37.vcf.gz
    curl -OkL https://ftp.ensembl.org/pub/current_variation/Enformer/enformer_grch37.vcf.gz.tbi
    curl -OkL https://ftp.ensembl.org/pub/current_variation/Enformer/enformer_grch38.vcf.gz
    curl -OkL https://ftp.ensembl.org/pub/current_variation/Enformer/enformer_grch38.vcf.gz.tbi
    mv enformer_grch37.vcf.gz ${params.output}
    mv enformer_grch37.vcf.gz.tbi ${params.output}
    mv enformer_grch38.vcf.gz ${params.output}
    mv enformer_grch38.vcf.gz.tbi ${params.output}
    """
}

process PRIMATEAI {
    publishDir "${params.output}", mode: 'move'

    input:
    path 'PrimateAI_scores_v0.2.tsv.gz'
    path 'PrimateAI_scores_v0.2_hg38.tsv.gz'

    output:
    path 'PrimateAI_scores_v0.2_GRCh37_sorted.tsv.gz'
    path 'PrimateAI_scores_v0.2_GRCh37_sorted.tsv.gz.tbi'
    path 'PrimateAI_scores_v0.2_GRCh38_sorted.tsv.gz'
    path 'PrimateAI_scores_v0.2_GRCh38_sorted.tsv.gz.tbi'

    script:
    """
    gunzip -cf PrimateAI_scores_v0.2.tsv.gz \\
      | sed '12s/.*/#&/' | sed '/^\$/d' \\
      | awk 'NR<12{print \$0;next}{print \$0 | "sort -k1,1 -k 2,2n -V"}' \\
      | bgzip > PrimateAI_scores_v0.2_GRCh37_sorted.tsv.gz
    tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh37_sorted.tsv.gz
    gunzip -cf PrimateAI_scores_v0.2_hg38.tsv.gz \\
      | sed '12s/.*/#&/' | sed '/^\$/d' \\
      | awk 'NR<12{print \$0;next}{print \$0 | "sort -k1,1 -k 2,2n -V"}' \\
      | bgzip > PrimateAI_scores_v0.2_GRCh38_sorted.tsv.gz
    tabix -s 1 -b 2 -e 2 PrimateAI_scores_v0.2_GRCh38_sorted.tsv.gz
    """
}

process UTRANNOTATER {
    publishDir "${params.output}", mode: 'move'

    output:
    path 'uORF_5UTR_GRCh37_PUBLIC.txt'
    path 'uORF_5UTR_GRCh38_PUBLIC.txt'
    
    script:
    """
    curl -OkL https://raw.githubusercontent.com/Ensembl/UTRannotator/master/uORF_5UTR_GRCh37_PUBLIC.txt
    curl -OkL https://raw.githubusercontent.com/Ensembl/UTRannotator/master/uORF_5UTR_GRCh38_PUBLIC.txt
    """
}

process LOEUF {
    publishDir "${params.output}", mode: 'move'

    output:
    path 'loeuf_dataset.tsv.gz'
    path 'loeuf_dataset.tsv.gz.tbi'

    script:
    """
    curl -OkL https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7334197/bin/41586_2020_2308_MOESM4_ESM.zip
    unzip -qq -o 41586_2020_2308_MOESM4_ESM.zip
    mv supplement/supplementary_dataset_11_full_constraint_metrics.tsv.gz ./
    gunzip -c supplementary_dataset_11_full_constraint_metrics.tsv.gz \\
      | (head -n 1 && tail -n +2  | sort -t\$'\\t' -k 76,76 -k 77,77n ) > loeuf_temp.tsv
    sed '1s/.*/#&/' loeuf_temp.tsv > loeuf_dataset.tsv
    bgzip loeuf_dataset.tsv
    tabix -f -s 76 -b 77 -e 78 loeuf_dataset.tsv.gz
    """
}


process LOUEF38 {
    container 'utsuno/images:crossmap'
    containerOptions '-v :/output'

    publishDir "${params.output}", mode: 'move'

    output:
    path 'loeuf_dataset_GRCh38.tsv.gz'
    path 'loeuf_dataset_GRCh38.tsv.gz.tbi'

    script:
    """
    cat 
    """
}

process DBNSFP {
    publishDir "${params.output}", mode: 'move'

    output:
    path 'dbNSFP${params.dbnsfp_version}_grch38.gz'
    path 


    script:
    """
    curl -OkL https://dbnsfp.s3.amazonaws.com/dbNSFP${params.dbnsfp_version}.zip
    unzip dbNSFP${params.dbnsfp_version}.zip

    zgrep -h -v ^#chr dbNSFP${params.dbnsfp_version}_variant.chr* \\
      | sort -k1,1 -k2,2n - \\
      | cat h - \\
      | bgzip -c > dbNSFP${params.dbnsfp_version}_grch38.gz
    tabix -s 1 -b 2 -e 2 dbNSFP${params.dbnsfp_version}_grch38.gz

    """
}


workflow {
    // ALPHAMISSENSE()
    // EMFORMER()
    // DOSAGESENSITIVITY()
    if (params.primateai_dir) {      
        primateai_grch37_file_ch = Channel.fromPath("${params.primateai_dir}/PrimateAI_scores_v0.2.tsv.gz")
        primateai_grch38_file_ch = Channel.fromPath("${params.primateai_dir}/PrimateAI_scores_v0.2_hg38.tsv.gz")
        PRIMATEAI(primateai_grch37_file_ch, primateai_grch38_file_ch)
    }
    // LOEUF()
    // REVEL()
    // UTRANNOTATER()
}