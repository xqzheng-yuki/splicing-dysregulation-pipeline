# Experimental Evaluation of Pipeline Modifications

In this document, we experimentally evaluate several modifications to the current pipeline to explore potential improvements. Our focus is on comparing different strategies for collecting read IDs with reduced false positives and selecting appropriate sequences for STAR junction identification.

## Task Tracker

- [ ] **(1)** Remove all lines marked as "fail to align" to reduce false positives.
- [ ] **(2)** Apply weighted alignment scoring to prioritize high-confidence reads.
- [ ] **(3)** Extract `m1` and `m2` read IDs to examine their patterns in gene-level calls.
- [ ] **(4)** Combine transcriptome-aligned reads from **primary alignments** with `d`-tagged reads, and remap using STAR to assess splicing site definitions.
- [ ] **(5)** Compare transcriptome-aligned reads from **primary alignments** vs. **all alignments** (no primary filter).
- [ ] **(6)** Based on (5), extract all transcriptome-aligned reads (no primary filter), combine with `d`-tagged reads, and remap with STAR to evaluate splicing site definitions.
- [ ] **(7)** Explore integration of `m1`/`m2` read data as additional input in the `d`-mapped pipeline.

```sh
samples=(CTX_104 CTX_108 CTX_128 CTX_154 CTX_120 CTX_125 CTX_147 CTX_148)
path="/mnt/gtklab01/xiaoqing/analysis/bam"
result_pre="/mnt/gtklab01/xiaoqing/analysis/list_gather"
salmon_pre="/mnt/gtklab01/xiaoqing/salmon"
```

```sh
for sample in "${samples[@]}"; do
    echo "Processed ${sample}" | tee -a ${result_pre}/report.txt
    # $base/analysis/bam/aligned/${sample}_d.bam -> lines that are actually aligned
    [ -e ${path}/aligned/${sample}_d.bam ] || samtools view -h ${path}/${sample}_d.bam | grep -v "AS:i:-2147483648" > ${path}/aligned/${sample}_d.bam
    # $base/analysis/list_gather/${sample}_mm_leftin.lst -> what unique in DECAP but align with mashmap region
    [ -e ${result_pre}/${sample}_mm_leftin.lst ] || samtools view ${path}/aligned/${sample}_d.bam | grep -v "intron" | cut -f1 > ${result_pre}/${sample}_mm_leftin.lst
    # append where are the non-overlap mashmap region appear in the pure mashmap aligned case
    echo "The composition of decoy mapped unique in decap mapped in mashmap:" >> ${result_pre}/report.txt
    grep -Ff ${result_pre}/${sample}_mm_leftin.lst ${salmon_pre}/${sample}/mashmap/aux_info/unmapped_names.txt | awk '$2 == "u" || $2 == "m1" || $2 == "m2" { count[$2]++ } END { for (k in count) print k, count[k] }' >> ${result_pre}/report.txt
    echo ""
    decap="${salmon_pre}/${sample}/full/aux_info/unmapped"
    mashmap="${salmon_pre}/${sample}/mashmap/aux_info/unmapped"
    [ -e ${decap}_names_sorted.txt ] || sort "${decap}_names.txt" > ${decap}_names_sorted.txt
    [ -e ${mashmap}_names_sorted.txt ] || sort "${mashmap}_names.txt" > ${mashmap}_names_sorted.txt
    echo ""
    echo "Tag | Overlap | DECAP_Count | mashmap_Count | %_in_DECAP | %_in_mashmap" | tee -a ${result_pre}/report.txt
    for tag in u m1 m2 d; do
        decap="${salmon_pre}/${sample}/full/aux_info/unmapped"
        mashmap="${salmon_pre}/${sample}/mashmap/aux_info/unmapped"
        [ -e ${decap}_${tag}.lst ] || grep $tag ${decap}_names.txt | cut -f1 -d' ' > ${decap}_${tag}.lst
        [ -e ${mashmap}_${tag}.lst ] || grep $tag ${mashmap}_names.txt | cut -f1 -d' ' > ${mashmap}_${tag}.lst
        
        count1=$(wc -l < ${decap}_${tag}.lst)
        count2=$(wc -l < ${mashmap}_${tag}.lst)

        overlap=$(comm -12 <(sort ${decap}_${tag}.lst) <(sort ${mashmap}_${tag}.lst) | wc -l)
        perc1=$(awk -v o=$overlap -v c=$count1 'BEGIN{printf("%.2f", o/c*100)}')
        perc2=$(awk -v o=$overlap -v c=$count2 'BEGIN{printf("%.2f", o/c*100)}')
        echo "$tag | $overlap | $count1 | $count2 | $perc1% | $perc2%" | tee -a ${result_pre}/report.txt
    done
    echo ""
    echo ">>> Non-overlap tag transitions (DECAP → mashmap):" >> ${result_pre}/report.txt
    for tag in u m1 m2 d; do
        echo "DECAP-$tag NOT in mashmap: appears in mashmap with tags:" >> ${result_pre}/report.txt
        decap="${salmon_pre}/${sample}/full/aux_info/unmapped"
        mashmap="${salmon_pre}/${sample}/mashmap/aux_info/unmapped"
        [ -e ${result_pre}/${sample}_decap_${tag}_unique.lst ] || comm -23 <(sort ${decap}_${tag}.lst) <(sort ${mashmap}_${tag}.lst) > ${result_pre}/${sample}_decap_${tag}_unique.lst
        # grep -Fxf "${result_pre}/${sample}_decap_${tag}_unique.lst" "${mashmap}_names.txt"
        grep -Ff ${result_pre}/${sample}_decap_${tag}_unique.lst "${mashmap}_names.txt" | awk '{count[$2]++} END {for (k in count) print k, count[k]}'  >> ${result_pre}/report.txt
        echo ""
    done
    echo ">>> Non-overlap tag transitions (mashmap → DECAP):" >> ${result_pre}/report.txt
    for tag in u m1 m2 d; do
        echo "mashmap-$tag NOT in DECAP: appears in DECAP with tags:" >> ${result_pre}/report.txt
        decap="${salmon_pre}/${sample}/full/aux_info/unmapped"
        mashmap="${salmon_pre}/${sample}/mashmap/aux_info/unmapped"
        [ -e ${result_pre}/${sample}_mashmap_${tag}_unique.lst ] || comm -23 <(sort ${mashmap}_${tag}.lst) <(sort ${decap}_${tag}.lst) > ${result_pre}/${sample}_mashmap_${tag}_unique.lst
        grep -Ff ${result_pre}/${sample}_mashmap_${tag}_unique.lst "${decap}_names.txt" | awk '{count[$2]++} END {for (k in count) print k, count[k]}' >> ${result_pre}/report.txt
        # awk 'NR==FNR{a[$1]; next} $1 in a' "${result_pre}/${sample}_mashmap_${tag}_unique.lst" "${decap}_names.txt" >> ${result_pre}/report.txt
    done
    echo "" | tee -a ${result_pre}/report.txt
done
```

```sh
for sample in "${samples[@]}"; do
    path="$salmon_pre/${sample}/full/aux_info"
    samtools view -h $path/primary_mapping.bam \
        | awk '$0 ~ /^@/ || $3 ~ /ENSMUST/' \
        | samtools view -bo $path/primary_trxn.inter.bam
    samtools sort --write-index -o $path/primary_trxn.bam##idx##$path/primary_trxn.bam.bai --output-fmt BAM $path/primary_trxn.inter.bam
    samtools view $path/primary_trxn.bam | cut -f1 > $path/trxn_id.lst
done
```

```sh

```
