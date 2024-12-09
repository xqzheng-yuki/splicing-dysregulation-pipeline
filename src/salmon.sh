#!/bin/bash
# set -o pipefail
# ## the actual name of the script directory
# SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# source "${SCRIPT_DIR}/atlas-config.sh"

# mkdir -p $SALMON_OUTPUT_DIR/$DATE/"mapping"

INDEX_DIR=$1
SAMPLE_NAME=$2
READS_1=$3
READS_2=$4
OUTPUT_DIR=$5

# activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon

# # read input read files from separate file
# echo $READS_FILE
# INDEX_DIR="/mnt/gtklab01/xiaoqing/salmon/index/decoy5/salmon_index"

# main body
# run salmon for each sample selective alignment (decoy mode)
# while read SAMPLE_NAME; do
    salmon quant -l A -i "${INDEX_DIR}" \
    -1 "${READS_1}" \
    -2 "${READS_2}" \
    --useVBOpt --writeUnmappedNames \
    --writeMappings="${OUTPUT_DIR}/${SAMPLE_NAME}/aux_info/mapping.bam" \
    -o "${OUTPUT_DIR}/${SAMPLE_NAME}" --threads 8
    # done < $READS_FILE

conda deactivate