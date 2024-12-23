#!/bin/bash

# File: separate_unmapped.sh
# Purpose: This script group sequences that tagged d, m1, m2, u into separate file.
# Dependencies: bash
# Input: unmapped_names.txt as the output of Salmon SA mode
# Output: Text files containing each type of unmapped sequence name, with each file saved as .lst
# Usage: bash separate_unmapped.sh
# Version History:
# v1.0 - Initial test release
# v2.0 - Initial general release
# Note: Ensure that the directory is the same or redirect by your need
set -o pipefail
## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "${SCRIPT_DIR}/atlas-config.sh"
PROJECT_DIR=$(dirname ${SCRIPT_DIR})

while read SAMPLE_NAME; do
    echo "start working on ${SAMPLE_NAME}"
    cd $SALMON_OUTPUT_DIR/$DATE/${SAMPLE_NAME}/"aux_info"
    INPUT_FILE="$SALMON_OUTPUT_DIR/$DATE/${SAMPLE_NAME}/aux_info/unmapped_names.txt"
    for UNMAPPED_TYPE in d m1 m2 u ; do
	    echo "working on $UNMAPPED_TYPE"
	    OUTPUT_FILE="./unmapped_${UNMAPPED_TYPE}.lst"
	    grep "${UNMAPPED_TYPE}$" $INPUT_FILE | cut -f1 -d" " > $OUTPUT_FILE
    done
done < $READS_FILE
