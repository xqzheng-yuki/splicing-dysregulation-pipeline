
output_dir = config['output_dir']
os.makedirs(output_dir,exist_ok=True)
container: "docker://continuumio/miniconda3"

rule intron_intervals:
    params:
        bedtools=config["bedtools"]
    input:
        genes=f"{output_dir}/genes.bed",
        exons=f"{output_dir}/exons.merged.bed",
    output:
        introns=f"{output_dir}/introns.bed",
    log:
        "logs/introns.log"
    shell:
        """
        {params.bedtools} subtract -s -a {input.genes} -b {input.exons} -nonamecheck > {output.introns}
        """

rule sort_introns:
    input:
        in_file=f"{output_dir}/introns.bed"
    output:
        f"{output_dir}/introns.sorted.bed"
    params:
        ## Add optional parameters for sorting order
    log:
        "logs/introns.sorted.log"
    wrapper:
        "v5.5.0/bio/bedtools/sort"

rule merge_introns:
    input:
        f"{output_dir}/introns.sorted.bed",
    output:
        f"{output_dir}/introns.merged.bed",
    params:
        extra="-s -c 4,5,6 -o distinct,distinct,distinct ",
    log:
        "logs/introns.merged.log",
    wrapper:
        "v5.5.0/bio/bedtools/merge"

rule exclude_exons:
    """
    Because all the operations are stranded, but alignment is not, let's exclude any
    intronic intervals that overlap with any exons.
    """
    input:
        left=f"{output_dir}/introns.merged.bed",
        right=f"{output_dir}/exons.merged.bed",
    output:
        f"{output_dir}/intronic.pure.bed",
    params:
        extra="-v"
    log:
        "logs/exclude_exons.log",
    wrapper:
        "v5.5.0/bio/bedtools/intersect"

rule intron_sequences:
    params:
        bedtools=config['bedtools']
    input:
        intervals=f"{output_dir}/intronic.pure.bed",
        fasta=get_genome_fasta,
    output:
        f"{output_dir}/intronic_found.fa",
    log:
        "logs/intron.seqs.log"
    shell:
        "{params.bedtools} getfasta -s -name -fi {input.fasta} -bed {input.intervals} -fo {output} 2> {log}"
    