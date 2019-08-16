def get_control(wildcards):
    return md.loc[(wildcards.sample),["Control"]]


rule bb2bed:
    input:
        "samples/bigBed/{sample}.all.bb"
    output:
        temp("samples/bed/{sample}.bed")
    run:
        bigBedToBed = config["bed_tool"]

        shell("""{bigBedToBed} {input} {output}""")

rule call_peaks:
    input:
        case = "samples/bed/{sample}.bed",
        control = get_control
    output:
        "results/motifs/{sample}/homerResults.html"
    params:
        assembly = config["assembly"],
        peak = get_peak
    conda:
        "../envs/peaks.yaml"
    shell:
        """scripts/callPeaks.sh -t {input.case} -c {input.control} -n {wildcards.sample} -p {params.peak} -a {params.assembly}"""

