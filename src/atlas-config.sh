#!/usr/bin/env bash

## top level
BASE_DIR="/mnt/gtklab01"
XQ_DIR="${BASE_DIR}/xiaoqing"
## This sets up the common configuration for files on atlas.

STAR_INPUT_DIR="${XQ_DIR}/star"

## filter fastq
RESULT_DIR="${XQ_DIR}/salmon/full_output"
STAR_OUTPUT_DIR="${XQ_DIR}/star/data"

## salmon
INDEX_DIR="${XQ_DIR}/salmon/decoy/index"
READS_DIR="${BASE_DIR}/linglab/tdp43/fastq"
SALMON_OUTPUT_DIR="${XQ_DIR}/salmon/full_output"
