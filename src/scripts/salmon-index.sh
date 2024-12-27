#!/bin/bash

# ## the actual name of the script directory
# SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# source "${SCRIPT_DIR}/atlas-config.sh"

# outfolder="${XQ_DIR}/salmon/index/decoy2"
# cd $outfolder

gentrome=$1
decoy=$2
output=$3

# activate conda environment
source ~e1124735/miniconda3/etc/profile.d/conda.sh
conda activate salmon
salmon index -t $gentrome --decoys $decoy -k 31 -p 8 -i $output
conda deactivate