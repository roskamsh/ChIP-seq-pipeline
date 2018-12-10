rule mapReads:
    input:
        reads = "samples/raw/{sample}.fastq.gz" if config["seq_type"]=="SE" else "samples/raw/{sample}_R1.fastq.gz|samples/raw/{sample}_R2.fastq.gz"
    output:
        "samples/bigBed/{sample}.mapped.reads.bb"
    params:
        assembly = config["assembly"],
        filter = config["filter"]
    conda:
        "../envs/bowtie.yaml"
    shell:
        """scripts/mapReads.sh -i {input.reads} -n {wildcards.sample} -a {params.assembly} -f {params.filter} -o samples/bigBed/"""

