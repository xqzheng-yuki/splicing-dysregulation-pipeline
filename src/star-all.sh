#!/bin/bash

# File: star-all.sh
# Purpose: This script align and map the unmapped read
# Dependencies: bash, STAR
# Input: STAR index file, compressed fastq file
# Output: STAR standard output files, with alignment file as "BAM file sorted by coordinate"
# Usage: sh star-all.sh
# Version History:
# v1.0 Initial release
# Note: Ensure that the directory is the same or redirect by your need

## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "${SCRIPT_DIR}/atlas-config.sh"

# directory
GENOME_DIR="/mnt/gtklab01/linglab/mmusculus_annotation_files/STAR_v2.7.9a_index_mmusculus_gencode.vM29"
READS_DIR="/mnt/gtklab01/xiaoqing/star/data"
OUTPUT_DIR="/mnt/gtklab01/xiaoqing/star"
# check if output directory exist
mkdir -p $OUTPUT_DIR

readarray -t samples < $READS_FILE 

ONEFILES=$(printf "%s_unmapped_seq_1.fq.gz," ${samples[@]})
ONEFILES=${ONEFILES::-1}
TWOFILES=$(printf "%s_unmapped_seq_2.fq.gz," ${samples[@]})
TWOFILES=${TWOFILES::-1}
RGLINE=$(printf "ID:%s , " ${samples[@]})
RGLINE=${RGLINE::-3}

cd $READS_DIR

STAR --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $ONEFILES $TWOFILES \
    --outSAMattrRGline $RGLINE \
    --outFileNamePrefix $OUTPUT_DIR/unmapped \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat --runThreadN 8 --outSAMattributes MD NH XS --outSAMunmapped Within --twopassMode Basic
