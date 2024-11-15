#!/bin/bash

## the actual name of the script directory
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
source "${SCRIPT_DIR}/atlas-config.sh"

outfolder="${XQ_DIR}/salmon/index/decoy2"
cd $outfolder

# activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate salmon
salmon index -t gentrome.fa --decoys decoysNintronic.txt -k 31 -p 8 -i salmon_index
