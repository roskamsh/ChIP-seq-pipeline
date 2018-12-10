rule mapReads:
    input:
        reads = "samples/raw/{sample}.fastq.gz" if config["seq_type"]=="SE" else "samples/raw/{sample}_R1.fastq.gz|samples/raw/{sample}_R2.fastq.gz"
    output:
        "samples/bigBed/{sample}.all.bb"
    params:
        assembly = config["assembly"],
        filter = config["filter"]
    conda:
        "../envs/chip.yaml"
    shell:
        """scripts/mapReads.sh -i {input.reads} -n {wildcards.sample} -a {params.assembly} -f {params.filter} -o samples/bigBed/"""

rule makeTracks:
    input:
        "samples/bigBed/{sample}.all.bb"
    output:
        "samples/bigwig/{sample}.bw"
    params:
        assembly = config["assembly"]
    conda:
        "../envs/chip.yaml"
    shell:
        """scripts/makeTracks.sh -i {input} -o {output} -g {params.assembly}"""
