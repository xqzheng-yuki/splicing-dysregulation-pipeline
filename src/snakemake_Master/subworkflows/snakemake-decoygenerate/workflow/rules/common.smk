from pathlib import Path

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
