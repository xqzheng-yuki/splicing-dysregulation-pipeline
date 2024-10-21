#!/bin/bash

# GOAL
# grab the actual sequence that unmapped


# directory
READS_DIR="/mnt/gtklab01/linglab/tdp43/fastq"
RESULT_DIR="/mnt/gtklab01/xiaoqing/salmon/full_output"
OUTPUT_DIR="/mnt/gtklab01/xiaoqing/star/data"
# check if output directory exist
mkdir -p $OUTPUT_DIR

# read input read files from separate file
READS_FILE="/mnt/gtklab01/xiaoqing/sample_name.txt"

AUX_DIR="aux_info"
NAME_SUF="unmapped_names.txt"
ID_SUF="unmapped_id.lst"

index=1
while read SAMPLE_NAME; do
    echo "Sample ${index}: Working on ${SAMPLE_NAME}."
    # modified file
    echo "Working on getting pure ID list."
    cut -d" " -f1 $RESULT_DIR/$SAMPLE_NAME/$AUX_DIR/$NAME_SUF > $RESULT_DIR/$SAMPLE_NAME/$AUX_DIR/$ID_SUF
    # actural extract
    echo "Working on extracting ${SAMPLE_NAME}_1."
    seqtk subseq $READS_DIR/${SAMPLE_NAME}_1.fq.gz $RESULT_DIR/$SAMPLE_NAME/$AUX_DIR/$ID_SUF > $OUTPUT_DIR/${SAMPLE_NAME}_unmapped_seq_1.fq
    echo "Finished extracting ${SAMPLE_NAME}_1."
    echo "Working on extracting ${SAMPLE_NAME}_2."
    seqtk subseq $READS_DIR/${SAMPLE_NAME}_2.fq.gz $RESULT_DIR/$SAMPLE_NAME/$AUX_DIR/$ID_SUF > $OUTPUT_DIR/${SAMPLE_NAME}_unmapped_seq_2.fq
    echo "Finished extracting ${SAMPLE_NAME}_2."
    echo ""
    index=$(expr $index + 1)
done < $READS_FILE