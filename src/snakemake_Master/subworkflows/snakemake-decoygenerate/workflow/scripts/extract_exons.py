#!/usr/bin/env python3

import sys
import gffutils

sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.output[0], "w")
transcript_type = snakemake.config["transcript_type"]

db = gffutils.FeatureDB(snakemake.input[0],keep_order=True)

for f in db.features_of_type('exon'):
    if (not transcript_type) or (transcript_type  == f['transcript_type'][0]):
        print(*[f.seqid,f.start,f.end,f['gene_id'][0],f.score,f.strand],sep='\t')

