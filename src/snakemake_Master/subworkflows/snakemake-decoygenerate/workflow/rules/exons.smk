output_dir = config['output_dir']

rule extract_exons:
    input:
        get_gtf_db,
    output:
        f"{output_dir}/exons.bed"
    conda:
        "../envs/gffutils.yaml"
    log:
        "logs/exons.log"
    script:
        "../scripts/extract_exons.py"

rule sort_exons:
    input:
        in_file=f"{output_dir}/exons.bed"
    output:
        f"{output_dir}/exons.sorted.bed"
    params:
        ## Add optional parameters for sorting order
    log:
        "logs/exons.sorted.log"
    wrapper:
        "v5.5.0/bio/bedtools/sort"

rule merge_exons:
    input:
        f"{output_dir}/exons.sorted.bed",
    output:
        f"{output_dir}/exons.merged.bed",
    params:
        extra="-s -c 4,5,6 -o distinct,distinct,distinct ",
    log:
        "logs/exons.merged.bed.log",
    wrapper:
        "v5.5.0/bio/bedtools/merge"