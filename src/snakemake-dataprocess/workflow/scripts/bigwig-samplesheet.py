#!/usr/bin/env python
# Generate a samplesheet containing both bigwig information and
# original sample information
# This can be used for later Shiny app

import sys
import pandas as pd
import itertools
import yaml
from snakemake.script import snakemake
sys.stderr = open(snakemake.log[0], "w")

## Get the sample info

sample_info = pd.read_csv(snakemake.config["samples"], sep="\t", dtype=str, comment="#").set_index("sample")

## Get the 
data = itertools.product(snakemake.params[0],snakemake.params[1])


df = pd.DataFrame(data,columns=['sample','tag'])
df = pd.DataFrame(data={'sample': df['sample'],
                   'tag': df['tag'],
                   "bw": snakemake.input.bigwig,
                   "bam": snakemake.input.bam,
                   },dtype=str).set_index("sample")

sample_info.join(df,on="sample").to_csv(snakemake.output[0],sep="\t")
