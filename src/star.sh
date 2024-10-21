#!/bin/bash

# directory
GENOME_DIR="/mnt/gtklab01/linglab/mmusculus_annotation_files/STAR_v2.7.9a_index_mmusculus_gencode.vM29"
READS_DIR="/mnt/gtklab01/xiaoqing/star/data"
OUTPUT_DIR="/mnt/gtklab01/xiaoqing/star"
# check if output directory exist
mkdir -p $OUTPUT_DIR

# read input read files from separate file
READS_FILE="/mnt/gtklab01/xiaoqing/sample_name.txt"

index=1
while read SAMPLE_NAME; do
    echo "Sample ${index}: Working on ${SAMPLE_NAME}."
    # main body
    STAR --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $READS_DIR/${SAMPLE_NAME}_unmapped_seq_1.fq.gz $READS_DIR/${SAMPLE_NAME}_unmapped_seq_2.fq.gz \
    --outFileNamePrefix $OUTPUT_DIR/$SAMPLE_NAME/${SAMPLE_NAME}_ \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat --runThreadN 8 --outSAMattributes MD NH XS --outSAMunmapped Within --twopassMode Basic
    echo ""
    index=$(expr $index + 1)
done < $READS_FILE
