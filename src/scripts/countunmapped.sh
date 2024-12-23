#!/bin/bash

## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "${SCRIPT_DIR}/atlas-config.sh"

NAME_SUF="aux_info/unmapped_names.txt"
while read SAMPLE_NAME; do
    echo "start working on ${SAMPLE_NAME}."
    while IFS= read -r line; do
        tag=$(echo "$line" | cut -d" " -f2)
        seqid=$(echo "$line" | cut -d" " -f1)
        echo "$seqid" >> "$SALMON_OUTPUT_DIR/$SAMPLE_NAME/"aux_info/unmapped"_${tag}.lst"
    done < $SALMON_OUTPUT_DIR/$SAMPLE_NAME/$NAME_SUF
done < $READS_FILE
