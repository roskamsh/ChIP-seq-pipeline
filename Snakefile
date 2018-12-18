__author__ = "Breeshey Roskams-Hieter"
__email__ = "roskamsh@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

import datetime
import sys
import os
import pandas as pd

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"config.yaml"
project_id = config["project_id"]
seq_type = config["seq_type"]
md = pd.read_table(config["samples"], index_col=["SampleID"], dtype=str)

if config["seq_type"]=="SE":
    CASES, = glob_wildcards("samples/cases/{sample}.fastq.gz")
else:
    CASES, = glob_wildcards("samples/cases/{sample}_R1.fastq.gz")

if config["seq_type"]=="SE":
    CONTROLS, = glob_wildcards("samples/controls/{sample}.fastq.gz")
else:
    CONTROLS, = glob_wildcards("samples/controls/{sample}_R1.fastq.gz")

## multiple samples may use the same control input/IgG files
CONTROLS_UNIQUE = list(set(CONTROLS))

SAMPLES = CASES + CONTROLS_UNIQUE

rule_dirs = ['mapReads','makeTracks','bb2bed','macs2']
for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
        log_out = os.path.join(os.getcwd(), 'logs', rule)
        os.makedirs(log_out)
        print(log_out)

def message(mes):
    sys.stderr.write("|--- " + mes + "\n")

for sample in SAMPLES:
    message("Sample " + sample + " will be processed")

rule all:
    input:
        expand("samples/bigBed/{sample}.all.bb", sample = SAMPLES),
        expand("samples/bigwig/{sample}.bw", sample = SAMPLES),
        expand("results/macs2/{sample}/{sample}_peaks.xls", sample = CASES)

include: "rules/align.smk"
include: "rules/peaks.smk"
