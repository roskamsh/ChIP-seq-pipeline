def get_gsize():
    if (config["assembly"] == "hg38") | (config["assembly"] == "hg19"):
        return 2.7e9
    elif (config["assembly"] == "mm9") | (config["assembly"] == "mm10"):
        return 1.87e9
    elif (config["assembly"] == "ce6") | (config["assembly"] == "ce10"):
        return 9e7
    elif (config["assembly"] == "dm3") | (config["assembly"] == "dm6"):
        return 1.2e8
    else:
        return "ERROR"

def get_control(wildcards):
    return md.loc[(wildcards.sample),["Control"]]


rule bb2bed:
    input:
        "samples/bigBed/{sample}.all.bb"
    output:
        temp("samples/bed/{sample}.bed")
    run:
        bigBedToBed = config["bed_tool"]

        shell("bigBedToBed {input} {output}")

rule call_peaks:
    input:
        case = "samples/bed/{sample}.bed",
        control = get_control
    output:
        "results/macs2/{sample}/{sample}_peaks.xls" if get_peak == "narrow" else "results/SICER/{sample}/{sample}-W200-G600-islands-summary"
    params:
        assembly = config["assembly"],
        peak = get_peak
    conda:
        "../envs/peaks.yaml"
    shell:
        """scripts/callPeaks.sh -t {input.case} -c {input.control} -n {wildcards.sample} -p {params.peak} -a {params.assembly}"""

