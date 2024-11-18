#!/bin/bash
set -o pipefail
## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

source "${SCRIPT_DIR}/atlas-config.sh"

mkdir -p $SALMON_OUTPUT_DIR/$DATE

# activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon

# read input read files from separate file
echo $READS_FILE
INDEX_DIR="/mnt/gtklab01/xiaoqing/salmon/index/decoy2/salmon_index"
# main body
# run salmon for each sample selective alignment (decoy mode)
while read SAMPLE_NAME; do
    salmon quant -l A -i $INDEX_DIR \
    -1 $RNA_READS_DIR/${SAMPLE_NAME}_1.fq.gz \
    -2 $RNA_READS_DIR/${SAMPLE_NAME}_2.fq.gz \
    --useVBOpt --writeUnmappedNames \
    -o $SALMON_OUTPUT_DIR/$DATE/$SAMPLE_NAME --threads 8
    done < $READS_FILE

