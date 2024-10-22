#!/usr/bin/env python

from Bio import SeqIO
import gzip


fileA = "/mnt/gtklab01/xiaoqing/star/group_data/CTX_120_d_unmapped_seq_1.fq"
fileB = "/mnt/gtklab01/xiaoqing/star/data/CTX_120_unmapped_seq_1.fq.gz"

A_ids = set([ record.id for record in list(SeqIO.parse(fileA,"fastq"))])
print(len(A_ids))

with gzip.open(fileB,"rt") as handle:
    B_ids = set([ record.id for record in list(SeqIO.parse(handle,"fastq"))])
print(len(B_ids))

print(f"Are the IDs in A a subset of those in B? {A_ids.issubset(B_ids)}")
