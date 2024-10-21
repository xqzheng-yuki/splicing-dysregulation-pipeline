#!/bin/bash

# GOAL
# grab the actual sequence that unmapped

## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "${SCRIPT_DIR}/atlas-config.sh"

# check if output directory exist
mkdir -p $UNMAPPED_SEQ_DIR

AUX_DIR="aux_info"
NAME_SUF="unmapped_names.txt"
ID_SUF="unmapped_id.lst"

index=1
while read SAMPLE_NAME; do
    echo "Sample ${index}: Working on ${SAMPLE_NAME}."
    # modified file
    echo "Working on getting pure ID list."
    cut -d" " -f1 $SALMON_OUTPUT_DIR/$SAMPLE_NAME/$AUX_DIR/$NAME_SUF > $SALMON_OUTPUT_DIR/$SAMPLE_NAME/$AUX_DIR/$ID_SUF
    # actural extract
    echo "Working on extracting ${SAMPLE_NAME}_1."
    seqtk subseq $RNA_READS_DIR/${SAMPLE_NAME}_1.fq.gz $SALMON_OUTPUT_DIR/$SAMPLE_NAME/$AUX_DIR/$ID_SUF > $UNMAPPED_SEQ_DIR/${SAMPLE_NAME}_unmapped_seq_1.fq
    echo "Working on extracting ${SAMPLE_NAME}_2."
    seqtk subseq $RNA_READS_DIR/${SAMPLE_NAME}_2.fq.gz $SALMON_OUTPUT_DIR/$SAMPLE_NAME/$AUX_DIR/$ID_SUF > $UNMAPPED_SEQ_DIR/${SAMPLE_NAME}_unmapped_seq_2.fq
    echo ""
    index=$(expr $index + 1)
done < $READS_FILE
