rule mapReads_paired_cases:
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

rule mapReads_paired_controls:
    input:
        R1 = "samples/controls/{sample}_R1.fastq.gz",
        R2 = "samples/controls/{sample}_R2.fastq.gz"
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

rule mapReads_single_cases:
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

rule mapReads_single_controls:
    input:
        R1 = "samples/controls/{sample}.fastq.gz"
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
        assembly = config["assembly"],
	ext = config["extension"]
    conda:
        "../envs/chip.yaml"
    shell:
        """scripts/makeTracks.sh -i {input} -o {output} -g {params.assembly} -e {params.ext}"""

