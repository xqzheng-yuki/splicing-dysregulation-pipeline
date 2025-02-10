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
    filter_bw="/mnt/gtklab01/xiaoqing/2025-01-14/filter/bw/${sample%}_d.bw"
    samtools view -N ${lst} -b ${bam} > ${filter_re}
    samtools index ${filter_re}
    bamCoverage -b ${filter_re} -o ${filter_bw}
done

# coverage
# samtools view -N /mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/CTX_104_d.lst -F 256  mapping.bam | cut -f3 | uniq -c | grep -v "," | awk '{print $2 "\t" $1}' > count_d.tsv
# samtools view -N /mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/CTX_104_d.lst -F 256  mapping.bam | cut -f3 | uniq -c | grep "," | awk '{print $2 "\t" $1}' > count_d_o.tsv
samples=("CTX_104" "CTX_108" "CTX_120" "CTX_125" "CTX_128" "CTX_147" "CTX_148" "CTX_154")
tags=("d" "m1" "m2")
mkdir -p "/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/count"
for sample in "${samples[@]}"; do
    # pre-filter for primary alignment
    ## beware: the bam file are different from the one use for bw
        bam="/mnt/gtklab01/xiaoqing/2025-01-14/salmon/${sample}/aux_info/mapping.bam"
        filter_bam="/mnt/gtklab01/xiaoqing/2025-01-14/salmon/${sample}/aux_info/primary_mapping.bam"
        samtools sort ${bam} | samtools view -F 256 -b -o ${filter_bam} --write-index
    # get count for different tag
    for tag in "${tags[@]}"; do
        lst="/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/${sample}_${tag}.lst"
        tsv_name="/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/count/count_${sample}_${tag}.tsv"
        samtools view -N ${lst} ${filter_bam} | cut -f3 | uniq -c | grep -v "," | awk '{print $2 "\t" $1}' > ${tsv_name}
    done
done
