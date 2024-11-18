#!/bin/bash

# File: star.sh
# Purpose: This script align and map the unmapped read
# Dependencies: bash, STAR
# Input: STAR index file, compressed fastq file
# Output: STAR standard output files, with alignment file as "BAM file sorted by coordinate"
# Usage: sh star-all.sh
# Version History:
# v1.0 Initial release
# v2.0 Make change to accept sorted input
# Note: Ensure that the directory is the same or redirect by your need
set -o pipefail
## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "${SCRIPT_DIR}/atlas-config.sh"

# check if output directory exist
mkdir -p $UNMAPPED_BAM_DIR/group/$DATE
READS_TAG_FILE="/mnt/gtklab01/xiaoqing/sample_name_tag.txt"
readarray -t samples < $READS_TAG_FILE 

ONEFILES=$(printf "%s_unmapped_seq_1.fq," ${samples[@]})
ONEFILES=${ONEFILES::-1}
TWOFILES=$(printf "%s_unmapped_seq_2.fq," ${samples[@]})
TWOFILES=${TWOFILES::-1}
RGLINE=$(printf "ID:%s , " ${samples[@]})
RGLINE=${RGLINE::-3}

cd $UNMAPPED_SEQ_DIR/$DATE

STAR --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $ONEFILES $TWOFILES \
    --outSAMattrRGline $RGLINE \
    --outFileNamePrefix $UNMAPPED_BAM_DIR/group/$DATE/unmapped \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand cat --runThreadN 8 --outSAMattributes MD NH XS --outSAMunmapped Within --twopassMode Basic

cd $UNMAPPED_BAM_DIR/group/$DATE
# Split the BAM file into separate files based on read groups or other tags, 
# naming the output files dynamically using the format string.
samtools split -f '%*_%!.%.' unmappedAligned.sortedByCoord.out.bam