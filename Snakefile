__author__ = "Breeshey Roskams-Hieter"
__email__ = "roskamsh@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline"""

import datetime
import sys
import os

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

configfile:"config.yaml"
project_id = config["project_id"]

SAMPLES, = glob_wildcards("samples/raw/{sample}_R1.fq")

rule_dirs = ['alignToBigBed']
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
        expand("samples/bigBed/{sample}.mapped.reads.bb", sample = SAMPLES)

include: "rules/align.smk"
