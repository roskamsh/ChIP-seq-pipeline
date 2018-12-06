rule mapReads:
    input:
        SE = "samples/raw/{sample}.fastq.gz",
        PE = "samples/raw/{sample}_R1.fastq.gz|samples/raw/{sample}_R2.fastq.gz" 
    output:
        "samples/bigBed/{sample}.bb"
    params:
        assembly = config["assembly"],
        filter = config["filter"],
        type = config["type"]
    conda:
        "../envs/bowtie.yaml"
    run:
        if params.type == "SE":
            shell("scripts/mapReads.sh -i={input.SE} -n={wildcards.sample} -a={params.assembly} -f={params.filter} -o=samples/bigBed/")
        else:
            shell("scripts/mapReads.sh -i={input.PE} -n={wildcards.sample} -a={params.assembly} -f={params.filter} -o=samples/bigBed/")
