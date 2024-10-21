# Salmon decoy sequence 

## Description

## Usage
### 1.`salmon.sh`
- align and map RNA-seq read using Salmon(https://github.com/COMBINE-lab/salmon) with selective alignment mode to get unmapped sequence
- before run script `salmon.sh` please make sure a conda environment named "salmon" that included salmon existed.
    - `conda create -n salmon salmon`
- output: 
    1. Quantification File
    2. Auxiliary Files - **`unmapped_names.txt`**
### 2.`filter_fastq.sh`
- adopt seqtk (package from https://github.com/lh3/seqtk) for extracting sequences with unmapped_names
- to install seqtk
    ```bash
    git clone https://github.com/lh3/seqtk.git;
    cd seqtk; make
    ```
- output:
    - one sequence name per line for each names (eg. `CTX_104_unmapped_seq_1.fq`, `CTX_104_unmapped_seq_2.fq`)
### 3.`star-allw.sh`
- used pass STAR index
    STAR v2.7.9.a index - generated using GRCm39.primary_assembly.genome.fa
    > STAR --runMode genomeGenerate --runThreadN 20 --genomeDir /mnt/gtklab01/linglab/mmusculus_annotation_files/STARv2.7.9a_index_mmusculus_gencode.vM27 --genomeFastaFiles /mnt/gtklab01/linglab/mmusculus_annotation_files/GRCm39.primary_assembly.genome.fa --sjdbGTFfile /mnt/gtklab01/linglab/mmusculus_annotation_files/gencode.vM27.primary_assembly.annotation.gtf
    Folder: STARv2.7.9a_index_mmusculus_gencode.vM27
    Annotatioin used:
    - GENCODE Primary assembly GRCm39 mouse FASTA: GRCm39.primary_assembly.genome.fa
> STAR --runMode alignReads \
    --genomeDir $GENOME_DIR \
    --readFilesIn $ONEFILES $TWOFILES \
    --outSAMattrRGline $RGLINE \
    --outFileNamePrefix $OUTPUT_DIR/unmapped \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesCommand zcat --runThreadN 8 --outSAMattributes MD NH XS --outSAMunmapped Within --twopassMode Basic
- output: unmappedAligned.sortedByCoord.out.bam
further manipulate: split the bam file by sample
> samtools split -f '%*_%!.%.' unmappedAligned.sortedByCoord.out.bam

## Project structure


## Contact
Author: Zheng Xiaoqing
Email: e1124735@u.nus.edu
Mentor: Greg Tucker-Kellogg
Email: dbsgtk@nus.edu.sg