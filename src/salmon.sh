#!/bin/bash

# directory
INDEX_DIR="/mnt/gtklab01/xiaoqing/salmon/decoy/index"
READS_DIR="/mnt/gtklab01/linglab/tdp43/fastq"
OUTPUT_DIR="/mnt/gtklab01/xiaoqing/salmon/full_output"
# check if output directory exist
mkdir -p $OUTPUT_DIR

# activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon

# read input read files from separate file
READS_FILE="/mnt/gtklab01/xiaoqing/sample_name.txt"
# READS_FILE="/mnt/gtklab01/xiaoqing/test_name.txt"

# main body
# run salmon for each sample selective alignment (decoy mode)
while read SAMPLE_NAME; do
    salmon quant -l A -i $INDEX_DIR \
    -1 $READS_DIR/${SAMPLE_NAME}_1.fq.gz \
    -2 $READS_DIR/${SAMPLE_NAME}_2.fq.gz \
    --useVBOpt --writeUnmappedNames \
    -o $OUTPUT_DIR/$SAMPLE_NAME --threads 8
    done < $READS_FILE

# deactivate conda environment
conda deactivate