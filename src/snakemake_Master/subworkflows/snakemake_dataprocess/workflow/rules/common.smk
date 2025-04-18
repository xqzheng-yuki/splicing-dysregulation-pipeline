from snakemake.utils import validate
import pandas as pd
import yaml
from pathlib import Path

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype=str, comment="#").set_index(
    "sample", drop=False)
samples.index.names = ["sample_id"]

def drop_unique_cols(df):
    singular_cols = df.nunique().loc[(df.nunique().values <= 1)].index
    return df.drop(singular_cols, axis=1)
 
 
samples = drop_unique_cols(samples)
validate(samples, schema="../schemas/samples.schema.yaml")
 
wildcard_constraints:
    sample="|".join(samples.index)


def get_gtf_file(wildcards):
    return config["Paths"]["gtffile"]

def get_gtf_db(wildcards):
    if "gtfdb" not in config["Paths"]:
        config["Paths"]["gtfdb"] = f"{config['Paths']['gtffile']}.db"
    return config["Paths"]["gtfdb"]

def get_genome_fasta(wildcards):
     return config['Paths']['genomefile']

def get_chrom_sizes(wildcards):
    return config['Paths']['chromsize']


get_gtf_db("")


def salmon_extra(wildcards):
    output_dir = config['output_dir']
    "--useVBOpt --writeUnmappedName --writeMappings={output_dir}/{{sample}}/aux_info/mapping.bam"