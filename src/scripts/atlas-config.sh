#!/usr/bin/env bash

## top level
BASE_DIR="/mnt/gtklab01"
XQ_DIR="${BASE_DIR}/xiaoqing"
READS_FILE="${XQ_DIR}/sample/sample_name.txt"
G_DIR="${BASE_DIR}/linglab/mmusculus_annotation_files"
DATE=$(date +%b_%d)
## This sets up the common configuration for files on atlas.

STAR_INPUT_DIR="${XQ_DIR}/star"

## salmon
INDEX_DIR="${XQ_DIR}/salmon/index/decoy1"
RNA_READS_DIR="${BASE_DIR}/linglab/tdp43/fastq"
SALMON_OUTPUT_DIR="${XQ_DIR}/salmon/full_output"

## filter fastq
UNMAPPED_SEQ_DIR="${XQ_DIR}/star/group_data"

## star-all
# pre-establish genome index
GENOME_DIR="${G_DIR}/STAR_v2.7.9a_index_mmusculus_gencode.vM29"
FASTAfiles="${G_DIR}/GRCm39.primary_assembly.genome.fa"
GTFfiles="${G_DIR}/gencode.vM29.primary_assembly.annotation.gtf"
UNMAPPED_BAM_DIR="${XQ_DIR}/star/results"