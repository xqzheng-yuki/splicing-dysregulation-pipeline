# Salmon decoy sequence 

## Description

## Usage
### 1. start with Salmon
#### a. Create Decoy and Build Salmon Index
    `generateDecoyTranscriptome.sh`

    `salmon-index.sh`

    `salmon.sh`
    - align and map RNA-seq read using [Salmon](https://github.com/COMBINE-lab/salmon) with selective alignment mode to get unmapped sequence
    - before run script `salmon.sh` please make sure a conda environment named "salmon" that included salmon existed.
        - `conda create -n salmon salmon`
    - input decoy index retrieved from [refgenie](http://refgenomes.databio.org/v3/genomes/splash/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1)

    | Asset name:tag  | Asset Description | Asset/achive size |
    | :---------------- | :------ | :---- |
    | salmon_partial_sa_index:default | Transcriptome index for salmon, produced with salmon index using partial selective alignment method. <br> Preparation includes transcriptome mapping to the genome and extraction of the relevant portion out from the genome and indexing it along with the transcriptome.<br>Recipe source -- https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh | 1.9GB / 802.1MB    |

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
### 3.`star-all.sh`
- Usage
    > sh star-all.sh
- Main Output: unmappedAligned.sortedByCoord.out.bam
- Notice:
    - used STAR index generate for previous project
    STAR v2.7.9.a index - generated using GENCODE Primary assembly GRCm39 mouse FASTA (GRCm39.primary_assembly.genome.fa)
    Folder: /mnt/gtklab01/linglab/mmusculus_annotation_files/STAR_v2.7.9a_index_mmusculus_gencode.vM29
- further manipulate: split the bam file by sample
    > samtools split -f '%*_%!.%.' unmappedAligned.sortedByCoord.out.bam

## Project structure


## Contact
Author: Zheng Xiaoqing
Email: e1124735@u.nus.edu
Mentor: Greg Tucker-Kellogg
Email: dbsgtk@nus.edu.sg
