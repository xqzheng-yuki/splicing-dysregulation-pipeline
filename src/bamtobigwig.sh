#!/bin/bash

## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source "${SCRIPT_DIR}/atlas-config.sh"

cd $UNMAPPED_BAM_DIR

for file in unmappedAligned.sortedByCoord.out_CTX*; do
    part_of_filename=$(echo "$file" | grep -o 'CTX_[^\.]*')
    echo "Processing file: $part_of_filename"
    samtools index $file
    output_name="unmapped"_$part_of_filename".bw"
    bamCoverage -b $file -o $output_name
done
