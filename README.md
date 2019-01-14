# ChIP-seq-pipeline
This is a pipeline to run basic Chromatin immunoprecipitation sequencing analysis for single-end or paired-end data.

This is a package of bash scripts that enable reading, processing and analysis of ChIP-seq Omics' datasets. This package implements the Snakemake management workflow system and is currently implemented to work with the cluster management and job scheduling system SLURM. This snakemake workflow utilizes conda installations to download and use packages for further analysis, so please ensure that you have installed miniconda prior to use.

# Questions/issues
Please add an issue to the cfRNA-seq repository. We would appreciate if your issue included sample code/files (as appropriate) so that we can reproduce your bug/issue.

# Contributing
We welcome contributors! For your pull requests, please include the following:

Sample code/file that reproducibly causes the bug/issue
Documented code providing fix
Unit tests evaluating added/modified methods.

# Usage:

Locate raw files.

```
$ cd /path/to/raw/data
$ ls -alh
```

Check md5sum.

```
$ md5sum –c md5sum.txt > md5sum_out.txt
```

Move your files into the archive to be stored.

```
$ mv /path/to/raw/data /path/to/archive
```

Check md5sum again in the archive to ensure your sequencing files are not corrupted.

```
$ md5sum –c md5sum.txt > md5sum_out.txt
```

Clone the ChIP-seq Pipeline into your working directory.

```
$ git clone https://github.com/ohsu-cedar-comp-hub/ChIP-seq-pipeline.git
```

Create a `samples/raw` directory, a `logs` directory and a `data` directory (if they do not exist) in your `wdir()`.

```
$ mkdir logs
$ mkdir data
$ mkdir samples
$ mkdir samples/raw
```

Symbollically link the fastq files of your samples to the `wdir/samples/raw` directory using a bash script loop in your terminal.

* For PE data:
```
$ cd samples/raw
$ ls -1 /path/to/data/LIB*R1*gz | while read gz; do
    R1=$( basename $gz | cut -d _ -f 2 | awk '{print $1"_R1.fastq.gz"}' )
    R2=$( basename $gz | cut -d _ -f 2 | awk '{print $1"_R2.fastq.gz"}' )
    echo $R1 : $R2
    ln -s $gz ./$R1
    ln -s ${gz%R1_001.fastq.gz}R2_001.fastq.gz ./$R2
done
```
* For SE data:
```
$ cd samples/raw
$ ls -1 /path/to/data/LIB*gz | while read gz; do 
  R=$( basename $gz | cut -d '_' -f 2 | awk '{print $1".fastq.gz"}' )
  echo $R
  ln -s $gz ./$R
done
```

Upload your metadata file to the `data` directory, with the correct columns:
* This file should be formatted as such:
```
SampleID  Control   Peak_call
H3K4Me3-1  samples/bed/IgG.bed  narrow
```
* Each row should be a `CASE` (i.e. not an Input/IgG), with the `Control` for that sample listed (written as `samples/bed/{ControlID}.bed`) and the type of peak calling you would like (either `broad` or `narrow`)
* All elements in this file should be tab-separated
* An example file is located in this repository here: `data/metadata.txt`

Edit the `omic_config.yaml` in your `wdir()`:
* Write the full path towards your `metadata.txt` file under the `samples` category
* Specify the genome `assembly` you would like to use
    * Current options for this are: hg19 and hg38
* Specify your `seq_type`
    * Either `SE` or `PE`

Do a dry-run of snakemake to ensure proper execution before submitting it to the cluster (in your wdir).

```
$ snakemake -np --verbose
```

Once your files are symbolically linked, you can submit the job to exacloud via your terminal window.

```
$ sbatch submit_snakemake.sh
```

To see how the job is running, look at your queue.

```
$ squeue -u your_username
```
