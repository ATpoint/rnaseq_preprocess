# rnaseq_preprocess

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000&logo=data%3Aimage%2Fjpeg%3Bbase64%2C%2F9j%2F4AAQSkZJRgABAQAAkACQAAD%2F2wBDABwcHBwcHDAcHDBEMDAwRFxEREREXHRcXFxcXHSMdHR0dHR0jIyMjIyMjIyoqKioqKjExMTExNzc3Nzc3Nzc3Nz%2F2wBDASIkJDg0OGA0NGDmnICc5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ubm5ub%2FwAARCACAAHgDASIAAhEBAxEB%2F8QAGgAAAgMBAQAAAAAAAAAAAAAABAUAAgMBBv%2FEADMQAAIBAwIDBgYBAwUAAAAAAAECAAMEESExEpHRBRVBUVJxEyIyYYGxQhQzwUNTgpLh%2F8QAGAEAAwEBAAAAAAAAAAAAAAAAAAIDAQT%2FxAAfEQEBAQACAgMBAQAAAAAAAAAAAQIDESExEkFRIjL%2F2gAMAwEAAhEDEQA%2FADe5rX1PzHSTua19T8x0jYbTsAUdzWvqfmOknc1r6n5jpG8yq1kpDLn8eMAW9zWvqfmOko3ZNmmrOw9yOk0qXdR9F%2BUfbeCHU5OspOP9SvJ%2BIbHs4f6jn2x0lk7OsnVmDVMLvt0lYZb%2FANqr7CG8yS0Z3bQw7OsD%2FNx746TVeybNvpdz%2BR0kk%2B4nLOSn7W7mtfU%2FMdJO5rX1PzHSbJcOujfMIalRagyplM6lb2WdzWvqfmOknc1r6n5jpG8kZpR3Na%2Bp%2BY6SRvJAODadnBtMq9UUU4vHwEILelLi4FEYGrGKWZnPExyTOMxYlmOSZtRoNWbTQDcy8kzHPbdVkqM54UGTDqdid6h%2FA6w6nTSmvCgxLxLu%2FSk459h1taA%2Fjn31mop01BCqADvLznEvmIndp%2BoyNCkf449pg9r4ofwYZkHYzsS4lHRQyspwwwZFYqeJTgxo6K4wwzF1Wk1I66jwMjrHXllg2jWFQYOjTeJwSpBG4jOlUFRc%2BPjKY334rZWskkko1wbRPc1fi1DjYaCMq7%2FDolhvsPzEsrxz7S5L9NKdM1XCL4x0iLTUIuwgtlT4U%2BId2%2FUNi7vd6bjPU7Ud0poXc4UakmILjtaq5K244F8zqeXhO9r1yai242UcR9ztEsMw9rZ69aocu7H8n%2FGJjgeIEk0p06lVxTpqWY%2BAjMZj5dV09tP1N0ubin9FRh%2Bc%2FvM2fs%2B8QZamSPsQYHAHFn2jd1KyUWw%2FEcZIwceJ0noGUOpVtjEHY9HiqPXP8Rwj3OpnoZPTYUuhRipl6L8DjyOhhNymU4xuP1AJzWfGsOZJlRbjpgyS8vZgt63yIvmcxbvp5w69yXT2MFpqTUUHzE6M%2BnPvzo7VQqhR4CWkkkHQ8l2ln%2BtqZ%2B2PbEBjzti3IIul20Vv8RJKy%2BC1yF2d0bSr8Th4gRgjx%2FEEkgHsaF7bXGlNhnyOh5TC9sEuV40wKg2Pn9jPKw%2B37QuaGnFxr5N13i%2FH8b29BYUDb2qIww27e5hkGtbqndoXp5GDgg%2BBhMStcYcSlT4xRtpHEU1Bh2H3MlyRlF2p0ZfIySlr9Te0kbHoRtVdkxjxEzFZiQDiXrj5QfKCx2mckqp4lB85aAUcIylKgBDaYPjEdz2QQS1qdPS3%2BDL9sViPh0VODniOPttNLPtRKgFO5PC%2B3F4HoY079sIKlKpRbhqqUP36zOe3dadVOFwGU%2BeonjrhaaV6iUjlAcCNL2yxjJJJGAuzuGt7hXB0JCsPsek9hPCHY4nuV%2Bke0TTYtF7XNQMQMYz5Q52CqWPgIonNy6666V4537H29V6hIbGkk5aD5S3mZI%2BP8%2BS79%2BBDrxIRF8ZDaB1k4WyNjHK1oNkcB8IRFqsVIYeEPVhUXIgHkr6r8a6dxsDwj2H%2FALmCR1V7GqDWjUDfZtDzHSBP2fepvTJ9iDKSwoMEgcIJA8gTicmzW9wv1UnH%2FEyvwq3%2B2%2F8A1MYM5IStpdP9NJvyMfuH0eyKznNdgg8hqekzsA7G3NzcKP4qQzH22H5nr5jRoUrdPh0hgfv3l6jimvEZPWvs0ga6fAFMeOpgMszF2LNuZtb0%2BN8nZZyW%2FLTon8wdSTgphZJpJOmTrw564NpV1Drgyw2nZoLmUqcGdRyhyIZUphxrv5wJlKHDQA5KiuNOUvFoJGomy12H1awAySYCuh3yJf4qFS2dBvANJIMbqmNsmYPdO2i%2FL%2B4l5Mw8xaMqVVpjXfyi2pUao2W5ShJJydZZEZzwqJDW7rwrnMy4ql2CruY1poKahRK0qS0h5k7may3HjrzUt67SSSSUI4Np2J%2B%2Bbb0vyHWd75tvS%2FIdYA3lWUMMMMxV3zbel%2BQ6yd823pfkOsAMegRqmswKlfqGJl3zbel%2BQ6znfNqf4PyHWAazVf7VT2gR7VsjvTbkOsnetngrwPg77dZlnhsWnVVm%2BkZmY7TsRtTbkOs075tRsj8h1kZw%2FtVvJ%2BCUtWOtQ4%2BwhqoqDCjAirvm29L8h1k75tvS%2FIdZXOJPSd1abyRR3zbel%2BQ6yd823pfkOsYpvJE%2FfNt6X5DrJAP%2F2Q%3D%3D)](https://sylabs.io/docs/)  

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
