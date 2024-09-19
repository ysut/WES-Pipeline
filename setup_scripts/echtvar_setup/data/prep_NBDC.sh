#!/bin/bash

set -eu

readonly URI="setup_workflow/ShellScripts/prep_GenomeAsia.sh"
readonly BASE_URI="https://browser.genomeasia100k.org/service/web/download_files"
readonly SUFFIX="annot.cont_withmaf.vcf.gz"
readonly MD5_FILE="md5sum.info.txt"
readonly DOWNLOAD_LOG="download.log"
readonly OUTPUT_PREFIX="GenomeAsia_snpindels_GRCh37"

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
  expected_checksum=$1
  filename=$2
  if [ -f "$filename" ]; then
    actual_checksum=$(md5sum "$filename" | awk '{print $1}')
    if [ "$expected_checksum" == "$actual_checksum" ]; then
      return 0
    fi
  fi
  return 1
}

# Download all VCF files
function download_vcfs() {
  download_file "${BASE_URI}/${MD5_FILE}"


  remove_old_uris_file
  touch "${URIS_FILE}"
  for i in {1..22}; do
    echo "${BASE_URI}/${i}.substitutions.${SUFFIX}" >> "${URIS_FILE}"
  done

  echo "${BASE_URI}/All.indels.${SUFFIX}" >> "${URIS_FILE}"
  # echo "${BASE_URI}/md5sum.info.txt" >> "${URIS_FILE}"

  # Download the MD5 file first
  download_file "${BASE_URI}/${MD5_FILE}"
  
  # Read the MD5 file and process each line
  tail -n +2 "${MD5_FILE}" | while IFS= read -r line; do
    expected_checksum=$(echo "${line}" | awk '{print $2}')
    filename=$(echo "${line}" | awk '{print $1}')
    if check_md5sum "${expected_checksum}" "${filename}"; then
      echo "File ${filename} already exists and matches the checksum. Skipping download."
    else
      download_file "${BASE_URI}/${filename}"
    fi
  done
  
  rm  "$URIS_FILE"
}

# Compare MD5 checksums
function compare_checksums() {
  if [ ! -f "$MD5_FILE" ]; then
    echo "MD5 checksum file not found! Now downloading ..."
    download_file "${BASE_URI}/${MD5_FILE}"
  fi

  tail -n +2 "${MD5_FILE}" | while IFS= read -r line; do
    expected_checksum=$(echo "${line}" | awk '{print $2}')
    filename=$(echo "${line}" | awk '{print $1}')
    if [ -f "${filename}" ]; then
      actual_checksum=$(md5sum "${filename}" | awk '{print $1}')
      if [ "${expected_checksum}" != "${actual_checksum}" ]; then
        echo "Checksum mismatch for ${filename}: expected ${expected_checksum}, got ${actual_checksum}" >> "${DOWNLOAD_LOG}"
      else
        echo "Checksum match for ${filename}"
      fi
    else
      echo "File ${filename} not found!" >> "${DOWNLOAD_LOG}"
    fi
  done
}

function preprocess() {
  if [ -s "${DOWNLOAD_LOG}" ]; then
    echo "Errors found in ${DOWNLOAD_LOG}. Exiting."
    exit 1
  fi
  echo "All files have been donwloaded."  

  echo "Decompressing ..."
  find -type f -name "*.vcf.gz" \
    | xargs -P 4 -I {} sh -c "gunzip {}"
  
  echo "Compressing as bgzip ..."
  find -type f -name "*.vcf" \
    | xargs -P 4 -I {} sh -c "bgzip -@ 2 {}"
  
  echo "Indexing with csi ..."
  find -type f -name "*.vcf.gz" \
    | xargs -P 4 -I {} sh -c "bcftools index {}" 
}

function concatenate() {
  echo ""
  bcftools concat \
    --allow-overlaps \
    --output "${OUTPUT_PREFIX}.vcf" \
    --output-type v \
    --threads 4 \
    *.substitutions.annot.cont_withmaf.vcf.gz All.indels.annot.cont_withmaf.vcf.gz
}

function compressvcf() {
  bcftools sort \
    --output "${OUTPUT_PREFIX}.sort.vcf" \
    "${OUTPUT_PREFIX}.vcf"

  bgzip -@ 4 "${OUTPUT_PREFIX}.sort.vcf"
  bcftools index "${OUTPUT_PREFIX}.sort.vcf.gz"
}

function main() {
  : > "${DOWNLOAD_LOG}"
  download_vcfs
  compare_checksums
  preprocess 
  concatenate
  compressvcf 
}

main
