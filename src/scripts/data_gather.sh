#!/bin/bash

bash lst_gather.sh

# comm -23 <(sort CTX_104-unmapped_d.lst) <(sort CTX_104_mm-unmapped_d.lst) > analysis/CTX_104_d.lst
for sample in CTX_*; do
    prefix="${sample%-unmapped_d.lst}"
    if [[ -f "${prefix}-unmapped_d.lst" && -f "${prefix}_mm-unmapped_d.lst" ]]; then
        comm -23 <(sort "${prefix}-unmapped_d.lst") <(sort "${prefix}_mm-unmapped_d.lst") > "analysis/${prefix}_d.lst"
        echo "Processed ${prefix}"
    else
        echo "Skipping ${prefix} due to missing files"
    fi
done

# samtools view -N /mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/CTX_104_d.lst d.sorted.bam
samples=("CTX_104" "CTX_108" "CTX_120" "CTX_125" "CTX_128" "CTX_147" "CTX_148" "CTX_154")
for sample in "${samples[@]}"; do
    bam="/mnt/gtklab01/xiaoqing/2025-01-14/star/result/${sample%}/d.sorted.bam"
    lst="/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/${sample%}_d.lst"
    filter_re="/mnt/gtklab01/xiaoqing/2025-01-14/filter/${sample%}_d.bam"
    samtools view -N ${lst} -b ${bam} > ${filter_re}
    samtools index ${filter_re}
done