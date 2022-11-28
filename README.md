# rnaseq_preprocess

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)  

## Introduction

**rnaseq_preprocess** is a Nextflow pipeline for RNA-seq quantification with `salmon`. The processing steps are `fastqc` first, then quantification with `salmon`, aggregation to gene level with `tximport` and a small summary report with `MultiQC`. Multiple fastq files per sample are supported. These technical replicates will be merged prior to quantification. Optional trimming to a fixed read length is possible. The pipeline is containerized via Docker and Singularity. Outputs can be found in `rnaseq_preprocess_results/` including command lines and software versions. The expected Nextflow version is 21.10.6.

Run the test profile to see which output is being produced. Downloading the Docker image may take a minute or two:

```bash
NXF_VER=21.10.6 nextflow run atpoint/rnaseq_preprocess -r main -profile docker,test_with_existing_idx,test_resources
```

## Details

**Indexing**

The indexing step must be run first and separately using the `--only_idx` flag. For this we need a reference transcriptome (gzipped), a reference genome as decoy (gzipped) and a GTF annotation file (gzipped).

`--only_idx`: trigger the indexing process  
`--idx_name`: name of the produced index, default `idx`  
`--idx_dir`: name of the directory inside `rnaseq_preprocess_results/` storing the index, default `salmon_idx`  
`--idx_additional`: additional arguments to `salmon index` beyond the defaults which are `--no-version-check -t -d -i -p --gencode`  
`--txtome`: path to the gzipped transcriptome fasta  
`--genome`: path to the gzipped genome fasta  
`--gtf`: path to the gzipped GTF file  
`--transcript_id`: name of GTF column storing transcript ID, default `transcript_id`  
`--transcript_name`: name of GTF column storing transcript name, default `transcript_name`  
`--gene_id`: name of GTF column storing gene ID, default `gene_id`  
`--gene_name`: name of GTF column storing gene name, default `gene_name`  
`--gene_type`: name of GTF column storing gene biotype, default `gene_type`  

For the indexing process, 30GB of RAM and 6 CPUs are required/hardcoded. On our HPC we use:  

```bash
NXF_VER=21.10.6 nextflow run atpoint/rnaseq_preprocess -r main  -profile singularity,slurm --only_idx \
    --genome path/to/genome.fa.gz --txtome path/to/txtome.fa.gz --gtf path/to/foo.gtf.gz \
    -with-report indexing_report.html -with-trace indexing_report.trace -bg > indexing_report.log
```    

**Quantification/tximport**

The pipeline runs via a [samplesheet](./test/samplesheet.csv) which is a CSV file with the columns:
`sample,r1,r2,libtype`. The first column is the name of the sample, followed by the paths to the R1 and
R2 files and the salmon [libtype](https://salmon.readthedocs.io/en/latest/library_type.html). If R2 is left blank
then single-end mode is triggered for that sample. Multiple fastq files (lane/technical replicates) are supported.
These must have the same sample column and will then be merged prior to quantification. Optionally, a `seqtk` module can
trim reads to a fixed read length, triggered by `--trim_reads` with a default of 75bp, controlled by `--trim_length`. 
The quantification then runs with the salmon options `--gcBias --seqBias --posBias` (for single-end without `--gcBias`). 
Transcript abundance estimates from `salmon` are then summarized to the gene level using [tximport](https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#Salmon) with its `lengthScaledTPM` option. That means returned gene-level counts are already corrected for average transcript length and can go into any downstream DEG analysis, for example with `limma`. Both a matrix of counts and effective gene lengths is returned.

Other options:

`--idx`: path to the salmon index folder  
`--tx2gene`: path to the tx2gene map matching transcripts to genes  
`--samplesheet`: path to the input samplesheet  
`--trim_reads`: logical, whether to trim reads to a fixed length  
`--trim_length`: numeric, length for trimming  
`--quant_additional`: additional options to `salmon quant` beyond `--gcBias --seqBias --posBias`  

We hardcoded 30GB RAM and 6 CPUs for the quantification. On our HPC we use:

```bash
NXF_VER=21.10.6 nextflow run atpoint/rnaseq_preprocess -r main -profile singularity,slurm \
    --idx path/to/idx --tx2gene path/to/tx2gene.txt --samplesheet path/to/samplesheet.csv \
    -with-report quant_report.html -with-trace quant_report.trace -bg > quant_report.log
```

**Other options**

`--merge_keep`: logical, whether to keep the merged fastq files  
`--merge_dir`: folder inside the output directory to store the merged fastq files  
`--trim_keep`: logical, whether to keep the trimmed fastq files  
`--trim_dir`: folder inside the output directory to store the trimmed fastq files  
`--skip_fastqc`: logical, whether to skip `fastqc`  
`--only_fastqc`: logical, whether to only run `fastqc` and skip quantification  
`--skip_multiqc`: logical, whether to skip `multiqc`  
`--skip_tximport`: logical, whether to skip the `tximport` process downstream of the quantification  
`--fastqc_dir`: folder inside the output directory to store the fastqc results  
`--multiqc_dir`: folder inside the output directory to store the multiqc results  
