#!/usr/bin/env bash

## top level
BASE_DIR="/mnt/gtklab01"
XQ_DIR="${BASE_DIR}/xiaoqing"
READS_FILE="${XQ_DIR}/sample_name.txt"
## This sets up the common configuration for files on atlas.

STAR_INPUT_DIR="${XQ_DIR}/star"

## salmon
INDEX_DIR="${XQ_DIR}/salmon/decoy/index"
RNA_READS_DIR="${BASE_DIR}/linglab/tdp43/fastq"
SALMON_OUTPUT_DIR="${XQ_DIR}/salmon/full_output"

## filter fastq
RESULT_DIR="${XQ_DIR}/salmon/full_output"
UNMAPPED_UNISEQ_DIR="${XQ_DIR}/star/data"
UNMAPPED_SEQ_DIR="${XQ_DIR}/star/group_data"

## star-all
# pre-establish genome index
GENOME_DIR="${BASE_DIR}/linglab/mmusculus_annotation_files/STAR_v2.7.9a_index_mmusculus_gencode.vM29"
UNMAPPED_BAM_DIR="${XQ_DIR}/star/results"