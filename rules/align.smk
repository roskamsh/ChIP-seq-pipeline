rule mapReads_paired:
    input:
        R1 = "samples/cases/{sample}_R1.fastq.gz",
        R2 = "samples/cases/{sample}_R2.fastq.gz"
    output:
        "samples/bigBed/{sample}.all.bb"
    params:
        assembly = config["assembly"],
        filter = config["filter"],
        type = config["seq_type"]
    conda:
        "../envs/chip.yaml"
    shell:
        """scripts/mapReads.sh -i {input.R1} -I {input.R2} -t {params.type} -n {wildcards.sample} -a {params.assembly} -f {params.filter} -o samples/bigBed/"""

rule mapReads_single:
    input:
        R1 = "samples/cases/{sample}.fastq.gz"
    output:
        "samples/bigBed/{sample}.all.bb"
    params:
        assembly = config["assembly"],
        filter = config["filter"],
        type = config["seq_type"]
    conda:
        "../envs/chip.yaml"
    shell:
        """scripts/mapReads.sh -i {input.R1} -t {params.type} -n {wildcards.sample} -a {params.assembly} -f {params.filter} -o samples/bigBed/"""


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

