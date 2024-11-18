#!/bin/bash


## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "${SCRIPT_DIR}/atlas-config.sh"

# check if output directory exist
mkdir -p $UNMAPPED_SEQ_DIR/$DATE

AUX_DIR="aux_info"
NAME_SUF="unmapped_names.txt"
ID_SUF="unmapped_id.lst"

index=1
while read SAMPLE_NAME; do
    echo "Sample ${index}: Working on ${SAMPLE_NAME}."
    # actural extract
    cd $SALMON_OUTPUT_DIR/$SAMPLE_NAME/$AUX_DIR
    for file in ./unmapped_*.lst
    do
        type=$(echo "${file}" | cut -d"_" -f2 | cut -d"." -f1)
        echo "working on $type type."
        seqtk subseq $RNA_READS_DIR/${SAMPLE_NAME}_1.fq.gz $file > $UNMAPPED_SEQ_DIR/$DATE/${SAMPLE_NAME}_${type}_unmapped_seq_1.fq
        seqtk subseq $RNA_READS_DIR/${SAMPLE_NAME}_2.fq.gz $file > $UNMAPPED_SEQ_DIR/$DATE/${SAMPLE_NAME}_${type}_unmapped_seq_2.fq
    done
    index=$(expr $index + 1)
    sleep 1
done < $READS_FILE
