#!/bin/bash

set -eu

readonly DATE="20231004"
readonly FASTA_RELEASE="release-112"

readonly BASE_URI="https://rgc-research.regeneron.com/me/downloads"
readonly FASTA_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
readonly FASTA_URI="https://ftp.ensembl.org/pub/${FASTA_RELEASE}/fasta/homo_sapiens/dna/${FASTA_FILE}"

readonly SCRIPT_DIR=$(cd $(dirname $0); pwd)
readonly PREFIX="rgc_me_variant_frequencies"
readonly URIS_FILE="./uris.txt"
readonly CONCATENATED_VCF38="${PREFIX}_all_${DATE}_GRCh38.vcf"
readonly CHAIN_FILE_DIR="share/chain_files"
readonly FASTA_FILE_DIR="share/reference_fasta"

# Remove old uris.txt file if it exists
function remove_old_uris_file() {
  if [ -e "$URIS_FILE" ]; then
    echo "Old URI file exists."
    rm -rf "$URIS_FILE"
    echo "Removed old uris.txt"
  fi
}

# Download function using curl
function download_with_curl() {
  echo "Downloading with curl: $1"
  curl -OkL -# "$1"
}

# Download function using wget
function download_with_wget() {
  echo "Downloading with wget: $1"
  wget -nv "$1"
}

# Determine which download function to use
function download_file() {
  if command -v curl > /dev/null; then
    download_with_curl "$1"
  elif command -v wget > /dev/null; then
    download_with_wget "$1"
  else
    echo "Neither curl nor wget is installed."
    exit 1
  fi
}

# Download all VCF files
function download_vcfs() {
  remove_old_uris_file
  touch "${URIS_FILE}"

  # Add URIs for each chromosome
  for i in {1..22}; do
    echo "${BASE_URI}/${DATE}/${PREFIX}_chr${i}_${DATE}.vcf.gz" >> "${URIS_FILE}"
    echo "${BASE_URI}/${DATE}/${PREFIX}_chr${i}_${DATE}.vcf.gz.tbi" >> "${URIS_FILE}"
  done
  for i in {X,Y}; do
    echo "${BASE_URI}/${DATE}/${PREFIX}_chr${i}_${DATE}.vcf.gz" >> "${URIS_FILE}"
    echo "${BASE_URI}/${DATE}/${PREFIX}_chr${i}_${DATE}.vcf.gz.tbi" >> "${URIS_FILE}"
  done

  # Read URIS_FILE, check file existence and download if necessary
  while IFS= read -r line; do
    filename=$(basename "${line}")
    if [ ! -f "${filename}" ]; then
      download_file "${line}"
    fi
  done < "${URIS_FILE}"
}

function download_chain_files() {
  chainfile_prefix="https://hgdownload.soe.ucsc.edu/goldenPath"
  hg38ToHg19="${chainfile_prefix}/hg38/liftOver/hg38ToHg19.over.chain.gz"
  hg19ToHg38="${chainfile_prefix}/hg19/liftOver/hg19ToHg38.over.chain.gz"
  
  mkdir -p "${SlsCRIPT_DIR}/${CHAIN_FILE_DIR}" && cd $_
  if [ ! -f "${SCRIPT_DIR}/${CHAIN_FILE_DIR}/hg38ToHg19.over.chain.gz" ]; then
    echo "Downloading chain file (hg38ToHg19) ..."
    download_file "${hg38ToHg19}"
  fi
  if [ ! -f "${SCRIPT_DIR}/${CHAIN_FILE_DIR}/hg19ToHg38.over.chain.gz" ]; then
    echo "Downloading chain file (hg19ToHg38) ..."
    download_file "${hg19ToHg38}"
  fi
}

function concatenate_pass_sort() {
  if [ ! -f "${CONCATENATED_VCF38}" ]; then
    echo "Concatenate VCFs using bcftools ..."
    bcftools concat \
      --allow-overlaps \
      --output "${CONCATENATED_VCF38}" \
      --output-type v \
      --threads 4 \
      *.vcf.gz
    
    echo "Filter the concatenated VCF ..."
    bcftools filter \
      --output "${CONCATENATED_VCF38%.*}.passed.vcf" \
      --include 'FILTER="PASS"' \
      "${CONCATENATED_VCF38}"

    echo "Sort the passed VCF ..." 
    bcftools sort \
      --output "${CONCATENATED_VCF38%.*}.passed.sort.vcf" \
      "${CONCATENATED_VCF38%.*}.passed.vcf"
    
    rm -f "${CONCATENATED_VCF38}"
  fi
}

function download_grch38_fasta() {
  if [ ! -f "${FASTA_FILE}" ]; then
    echo "Downloading GRCh38 fasta file ..."
    download_file "${FASTA_URI}"
    gunzip "${FASTA_FILE}"
    bgzip -@4 "${FASTA_FILE%.*}"
  fi
}

# Final step
function sort_compress_index() {
  echo "Compressing and indexing VCFs (GRCh37) ..."
  bcftools sort \
    --output "${SCRIPT_DIR}/RGC-ME/GRCh37/${PREFIX}_all_${DATE}_GRCh37_remap.passed.sort.vcf" \
    "${SCRIPT_DIR}/RGC-ME/GRCh37/${PREFIX}_all_${DATE}_GRCh37_remap.passed.vcf"
  rm -f "${SCRIPT_DIR}/RGC-ME/GRCh37/${PREFIX}_all_${DATE}_GRCh37_remap.passed.vcf"
  bgzip -@ 4 ${SCRIPT_DIR}/RGC-ME/GRCh37/${PREFIX}_all_${DATE}_GRCh37_remap.passed.sort.vcf
  tabix -p vcf ${SCRIPT_DIR}/RGC-ME/GRCh37/${PREFIX}_all_${DATE}_GRCh37_remap.passed.sort.vcf.gz
  
  echo "Compressing and indexing VCFs (GRCh38) ..."
  bgzip -@ 4 ${SCRIPT_DIR}/RGC-ME/GRCh38/${CONCATENATED_VCF38%.*}.passed.sort.vcf
  tabix -p vcf ${SCRIPT_DIR}/RGC-ME/GRCh38/${CONCATENATED_VCF38%.*}.passed.sort.vcf.gz
}

function main() {
  # Download chain files (38 -> 37 and 37-> 38)
  mkdir -p "${SCRIPT_DIR}/${CHAIN_FILE_DIR}" && cd $_
  download_chain_files
  cd ${SCRIPT_DIR}

  # Download fasta (GRCh38)
  mkdir -p "${SCRIPT_DIR}/${FASTA_FILE_DIR}" && cd $_
  download_grch38_fasta
  cd ${SCRIPT_DIR}

  # Download RGC-ME VCFs and concatenate them with sorting
  mkdir -p ${SCRIPT_DIR}/RGC-ME/Original && cd $_
  download_vcfs
  concatenate_pass_sort
  mkdir -p ${SCRIPT_DIR}/RGC-ME/GRCh38
  mv ${SCRIPT_DIR}/RGC-ME/Original/${CONCATENATED_VCF38%.*}.passed.sort.vcf \
    ${SCRIPT_DIR}/RGC-ME/GRCh38
  cd ${SCRIPT_DIR}  
  
  # Remapping 38 to 37
  mkdir -p ${SCRIPT_DIR}/RGC-ME/GRCh37
  time docker run --rm \
    -v ${SCRIPT_DIR}/RGC-ME:/RGC-ME \
    -v ${SCRIPT_DIR}/${CHAIN_FILE_DIR}:/chain_files \
    -v ${SCRIPT_DIR}/${FASTA_FILE_DIR}:/reference_fasta \
    utsuno/images:crossmap \
    CrossMap vcf \
      /chain_files/hg38ToHg19.over.chain.gz \
      /RGC-ME/GRCh38/${CONCATENATED_VCF38%.*}.passed.sort.vcf \
      /reference_fasta/${FASTA_FILE} \
      /RGC-ME/GRCh37/${PREFIX}_all_${DATE}_GRCh37_remap.passed.vcf \
      --chromid a

  # Compress and index VCFs GRCh38 and GRCh37
  sort_compress_index

  rm -rf ${SCRIPT_DIR}/RGC-ME/Original/${CONCATENATED_VCF38%.*}.passed.sort.vcf
}

main
