#!/bin/bash

set -eu

readonly RELEASE="b156"
readonly BASE_URI="https://ftp.ncbi.nih.gov/snp/archive/${RELEASE}/VCF"

readonly grch37_vcf="${BASE_URI}/GCF_000001405.25.gz"
readonly grch37_tbi="${grch37_vcf}.tbi"
readonly grch37_vcf_md5="${grch37_vcf}.md5"
readonly grch37_tbi_md5="${grch37_tbi}.md5"
readonly grch38_vcf="${BASE_URI}/GCF_000001405.40.gz"
readonly grch38_tbi="${grch38_vcf}.tbi"
readonly grch38_vcf_md5="${grch38_vcf}.md5"
readonly grch38_tbi_md5="${grch38_tbi}.md5"

readonly basename37="$(basename ${grch37_vcf})"
readonly basename38="$(basename ${grch38_vcf})"

readonly SCRIPT_DIR=$(cd $(dirname $0); pwd)
readonly final_vcf37="${basename37%.*}.dbSNP-${RELEASE}.singleALT.common.chrmod.sort.vcf.gz"
readonly final_vcf38="${basename38%.*}.dbSNP-${RELEASE}.singleALT.common.chrmod.sort.vcf.gz"


# Download function using curl
function download_with_curl() {
  echo "Downloading with curl: $1"
  curl -Ok -# "$1"
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
# Check MD5 checksum of existing file
function check_md5sum() {
  if md5sum --check "$1" | grep -Ev 'OK$' > /dev/null; then
    echo "[ERROR] MD5 checksum failed. ($1)"
    exit 1
  fi
}

function download_rule_files() {
  # mkdir -p "${SCRIPT_DIR}/share/chrmod_rules" 
  if [ ! -f "chr_mod_rules_37.txt" ]; then
    # cd "${SCRIPT_DIR}/share/chrmod_rules"
    download_file https://raw.githubusercontent.com/ysut/WES-Pipeline/develop/resources/chr_mod_rules_37.txt
  fi
  if [ ! -f "chr_mod_rules_38.txt" ]; then
    # cd "${SCRIPT_DIR}/share/chrmod_rules"
    download_file https://raw.githubusercontent.com/ysut/WES-Pipeline/develop/resources/chr_mod_rules_38.txt
  fi
  # cd "${SCRIPT_DIR}"
}

# Download all VCF files
function download_vcfs() {
  # For GRCh37
  if [ ! -f "${basename37}" ]; then
    download_file "${grch37_vcf}"
    download_file "${grch37_tbi}"
  fi
  download_file "${grch37_vcf_md5}"
  download_file "${grch37_tbi_md5}"
  check_md5sum "$(basename ${grch37_vcf_md5})"
  check_md5sum "$(basename ${grch37_tbi_md5})"

  # For GRCh38
  if [ ! -f "${basename38}" ]; then
    download_file "${grch38_vcf}"
    download_file "${grch38_tbi}"
  fi  
  download_file "${grch38_vcf_md5}"
  download_file "${grch38_tbi_md5}"
  check_md5sum "$(basename ${grch38_vcf_md5})"
  check_md5sum "$(basename ${grch38_tbi_md5})"
}

function modify_header() {
  echo "Modifying header ..."
  bcftools view --header-only "${basename37}" --output header37.txt
  bcftools view --header-only "${basename38}" --output header38.txt
  sed -i 's/ID=RS,Number=1,Type=Integer/ID=RS,Number=1,Type=String/' header37.txt
  sed -i 's/ID=RS,Number=1,Type=Integer/ID=RS,Number=1,Type=String/' header38.txt
  
  bcftools reheader --threads 4 \
    --header header37.txt \
    --output "${basename37%.*}.rehead.vcf.gz" \
    "${basename37}"

  bcftools reheader --threads 4 \
    --header header38.txt \
    --output "${basename38%.*}.rehead.vcf.gz" \
    "${basename38}"
  
  rm -rf header37.txt header38.txt
}

function extract_edit_vcf() {
  if [ ! -f "${final_vcf37}" ]; then
    # echo "Extracting single ALT and COMMON variants ..."
    bcftools view --threads 4 --output-type u  \
      --include 'COUNT(ALT)=1 && INFO/COMMON=1' \
      "${basename37%.*}.rehead.vcf.gz" \
    | bcftools annotate --threads 4 --output-type u  \
      --rename-chrs chr_mod_rules_37.txt \
    | bcftools sort --max-mem 8G --output-type z \
      --output "${final_vcf37}" 
    tabix -p vcf "${final_vcf37}"
    rm -rf "${basename37%.*}.singleALT.common.chrmod.vcf"
    rm -rf "${basename37%.*}.rehead.vcf.gz"
  fi
  if [ ! -f "${final_vcf38}" ]; then
    bcftools view --threads 4 --output-type u \
      --include 'COUNT(ALT)=1 && INFO/COMMON=1' \
      "${basename38%.*}.rehead.vcf.gz" \
    | bcftools annotate --threads 4 --output-type u \
      --rename-chrs chr_mod_rules_38.txt \
    | bcftools sort --max-mem 8G --output-type z \
      --output "${final_vcf38}" 
    tabix -p vcf "${final_vcf38}"
    rm -rf "${basename38%.*}.singleALT.common.chrmod.vcf"
    rm -rf "${basename38%.*}.rehead.vcf.gz"
  fi
}

function main() {
  # prepare rule files
  download_rule_files
  
  # Download and process dbSNP VCFs
  download_vcfs
  modify_header
  extract_edit_vcf
}

main