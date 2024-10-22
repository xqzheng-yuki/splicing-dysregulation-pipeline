#!/usr/bin/env bash

## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "${SCRIPT_DIR}/atlas-config.sh"
PROJECT_DIR=$(dirname ${SCRIPT_DIR})

echo $PROJECT_DIR

SALMON_OUTPUT_DIR=${PROJECT_DIR}/test_results
mkdir -p $SALMON_OUTPUT_DIR

READS_FILE=${PROJECT_DIR}/test_data/reads_file.txt

while read SAMPLE_NAME; do
    echo "start working on ${SAMPLE_NAME}"
    OUTPUT_DIR="$SALMON_OUTPUT_DIR/$SAMPLE_NAME/aux_info"
    mkdir -p $OUTPUT_DIR
    INPUT_FILE="${PROJECT_DIR}/test_data/${SAMPLE_NAME}_part_unmapped_name.txt"
    ls $INPUT_FILE
    for UNMAPPED_TYPE in d m1 m2 u ; do
	echo "working on $UNMAPPED_TYPE"
	OUTPUT_FILE="$OUTPUT_DIR/unmapped_${UNMAPPED_TYPE}.lst"
	grep "${UNMAPPED_TYPE}$" $INPUT_FILE | cut -f1 -d" " > $OUTPUT_FILE
    done
done < $READS_FILE
