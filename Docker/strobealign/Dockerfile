FROM continuumio/miniconda3:main

RUN set -ex \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict \
    && conda create --name strobealign python=3.13 strobealign samtools --yes

RUN set -ex \
    && apt-get update > /dev/null \
    && apt-get install --no-install-recommends --yes \
        tini \
        > /dev/null \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["tini", "--"]
CMD [ "/bin/bash" ]


# Using "Homo_sapiens_assembly38.fasta" and "Homo_sapiens_assembly38.fasta.fai"
# The refernce file was downloaded from 
# https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0/

# wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
# wget -c https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai