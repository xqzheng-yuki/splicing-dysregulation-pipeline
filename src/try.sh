#!/bin/bash

READS_FILE="/mnt/gtklab01/xiaoqing/sample_name.txt"
RESULT_DIR="/mnt/gtklab01/xiaoqing/salmon/full_output"
NAME_SUF="aux_info/unmapped_names.txt"
while read SAMPLE_NAME; do
    echo ${SAMPLE_NAME}
    cut -d" " -f 2 $RESULT_DIR/$SAMPLE_NAME/$NAME_SUF | sort | uniq -c
done < $READS_FILE