#!/bin/bash

DORADO_VERSION="0.7.3"
DORADO="/usr/local/bin/dorado"
DNA_MODEL="/usr/local/share/dorado/models/"
output_dir="results"

${DORADO} basecall \
  --emit-fastq \
  > ${output_dir}