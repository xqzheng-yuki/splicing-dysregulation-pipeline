#!/bin/bash

bash lst_gather.sh
samples=("CTX_104" "CTX_108" "CTX_120" "CTX_125" "CTX_128" "CTX_147" "CTX_148" "CTX_154")
tags=("d" "m1" "m2")
mkdir -p "/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/count"
mkdir -p "/mnt/gtklab01/xiaoqing/2025-01-14/analysis"

# comm -23 <(sort CTX_104-unmapped_d.lst) <(sort CTX_104_mm-unmapped_d.lst) > analysis/CTX_104_d.lst
for sample in "${samples[@]}"; do
    echo "Processed ${sample}"
    for tag in "${tags[@]}"; do
        if [[ -f "${sample}-unmapped_${tag}.lst" && -f "${sample}_mm-unmapped_${tag}.lst" ]]; then
            comm -23 <(sort "${sample}-unmapped_${tag}.lst") <(sort "${sample}_mm-unmapped_${tag}.lst") > "analysis/${sample}_${tag}.lst"
        else
            echo "Skipping ${sample} due to missing files"
        fi
    done
done

# coverage
# samtools view -N /mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/CTX_104_d.lst -F 256  mapping.bam | cut -f3 | uniq -c | grep -v "," | awk '{print $2 "\t" $1}' > count_d.tsv
# samtools view -N /mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/CTX_104_d.lst -F 256  mapping.bam | cut -f3 | uniq -c | grep "," | awk '{print $2 "\t" $1}' > count_d_o.tsv

for sample in "${samples[@]}"; do
    # pre-filter for primary alignment
    ## beware: the bam file are different from the one use for bw
    bam="/mnt/gtklab01/xiaoqing/2025-01-14/salmon/${sample}/aux_info/mapping.bam"
    filter_bam="/mnt/gtklab01/xiaoqing/2025-01-14/salmon/${sample}/aux_info/primary_mapping.bam"
    if [ ! -e ${filter_bam} ]; then
        echo "${filter_bam} does not exist."
        samtools sort ${bam} | samtools view -F 256 -b -o ${filter_bam} --write-index
    fi
    # get count for different tag
    for tag in "${tags[@]}"; do
        lst="/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/${sample}_${tag}.lst"
        tsv_name="/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/count/count_${sample}_${tag}.tsv"
        samtools view -N ${lst} ${filter_bam} | cut -f3 | uniq -c | grep -v "," | awk '{print $2 "\t" $1}' >| ${tsv_name}
    done
done

# samtools view -N /mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/CTX_104_d.lst primary_mapping.bam
# samtools view -N $lst ../salmon/CTX_120/aux_info/primary_mapping.bam
mkdir -p "/mnt/gtklab01/xiaoqing/2025-01-14/analysis/bam"
mkdir -p "/mnt/gtklab01/xiaoqing/2025-01-14/analysis/bw"
for sample in "${samples[@]}"; do
    bam="/mnt/gtklab01/xiaoqing/2025-01-14/salmon/${sample}/aux_info/primary_mapping.bam"
    for tag in "${tags[@]}"; do
        lst="/mnt/gtklab01/xiaoqing/2025-01-14-list_filter/analysis/${sample%}_${tag}.lst"
        filter_re="/mnt/gtklab01/xiaoqing/2025-01-14/analysis/bam/${sample%}_${tag}.bam"
        filter_bw="/mnt/gtklab01/xiaoqing/2025-01-14/analysis/bw/${sample%}_${tag}.bw"
        samtools view -N ${lst} ${bam} -b -o ${filter_re} --write-index
        bamCoverage -b ${filter_re} -o ${filter_bw}
    done
done