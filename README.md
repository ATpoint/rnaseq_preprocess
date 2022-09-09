# rnaseq_preprocess

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>

## Introduction

**rnaseq_preprocess** is an automated preprocessing pipeline for RNA-seq data implemented using [Nextflow](https://www.nextflow.io/) which is fully containerized to take care of all required software and ensure reproducibility. It performs an initial QC with `fastqc`, then builds a genome-decoyed index against a reference transcriptome and quantifies the fastq files against it with [salmon](https://salmon.readthedocs.io/en/latest/salmon.html). Eventually, the individual quantifications are summarized to the gene level with `tximport` returing a single matrix of raw counts. The `fastqc` and mapping statistics are summarized with `MultiQC`. The pipeline supports automated merging of technical/lane/sequencing replicates and auto-detects whether data are single- or paired-end. The starting point for the pipeline is a samplesheet.

# Details

**Indexing**<br>

The indexing requires a transcriptome and genome in fasta format and a reference annotation in GTF format. All files must be gzipped and end with `.gz`. The flags to provide these files are `--genome`, `--txtome` and `--gtf`. Also, for creating a transcripts to gene (tx2gene) map for gene level summarization after mapping the pipeline needs `--transcript_id`, `--transcript_name`, `--gene_id`, `--gene_name` and `--gene_type` which take the colnames of the GTF that store threse attributes. The defaults assume GENCODE reference annotations, and if using them these options do not need to be provided. For this we also pass the `--gencode` flag to the `salmon` indexing process via `--idx_additional` like `--idx_additional '--gencode'`. Any other salmon indexing options can be passed via this flag as well such as non-standard kmer length, e.g. `idx_additional '--gencode -k 17'` which makes sense when quantifying very short reads such as the 2x25bp data from the ImmGen consortium as the default kmer length is 31, and it must be shorter than the read length. THe putput of the indexing will be in the folder:
`rnaseq_preprocess_results/salmonIdx/name` where "name" is the value passed via `--idx_name` which defaults to simply "idx". One can use the pipeline to only build the index and then exit. It makes sense to build the index first and then store it somewhere for re-use, e.g.:

```bash
NXF_VER=21.10.6 \
    nextflow run main.nf \
    --publishmode 'copy' \
    --only_idx \ # this makes the pipeline stop after indexing
    --genome path/to/genome.fa.gz --txtome path/to/txtome.fa.gz --gtf path/to/foo.gtf.gz \
    -with-report indexing_report.html \
    -with-trace indexing_report.log
```

After this finishes successfully one can delete the `work` directory as the results have been copied to the output folder. See the [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) on choices for `--publishmode` and what the pro/cons are per choice.

**QC/Quantification**

The quantification requires a samplesheet with four columns such as the one in the test folder:

```bash
sample,r1,r2,libtype
sample1,$baseDir/test/sample1_1.fastq.gz,$baseDir/test/sample1_2.fastq.gz,A
sample2,$baseDir/test/sample2_1.fastq.gz,$baseDir/test/sample2_2.fastq.gz,A
sample2,$baseDir/test/sample2a_1.fastq.gz,$baseDir/test/sample2a_2.fastq.gz,A
sample3,$baseDir/test/sample3.fastq.gz,,A
```

The first row is the header which must be `sample,r1,r2,libtype`. These four fields are:
- `sample`: The sample name, this can be any user-defined name (no whitespaces here!)    
- `r1/r2`: The paths to the fastq files belonging to this sample. If this is single-end then leave `r2` empty, the pipeline will then detect this and trigger single-end mode.
- `libtype`: The library type compatible with the `--libtype` argument of `salmon`. Default is automatic detection (`A`). For Illumina stranded paired-end data that is mostly `ISR`. See the [salmon docs](https://salmon.readthedocs.io/en/latest/library_type.html) for appropriate choices of your libraries. 

Technical/lane replicates can provided here as well by using the **same** sample name (column1) for the fastq files being a technical replicate. In the above example `sample2` has two fastq file pairs that are a technical replicate. These will be merged internally prior to quantification. As of version v2.3 on the fastq files of technical replicates can have the same name (in different input folders). In prior versions this caused a name collision error.

Custom command line arguments can be passed to `salmon quant` (which runs the quantification) by `--quant_additional`, e.g. `--quant_additional '--numGibbsSamples 8'` to produce inferential replicates. By default we run paired-end data with the `salmon` flags `--posBias --seqBias --gcBias` and single-end data with the same excluding `--gcBias` via that `--quant_additional` flag. Use `--quant_additional ''` to turn that off. In case of **end-tagged** data such as QuantSeq one should use `--quant_additional '--noLengthCorrection'`.

Prior to quantification a `fastqc` check is run. If only this shall be executed without indexing and quantification, e.g. to check whether adapter contamination is present that require trimming then use `--only_fastqc`.

Once the quantification is done the pipeline will produce a gene-level count matrix for all samples using `tximport`. This count matrix can directly be used with applications such as `DESeq2` and `limma`. It uses the `lengthScaledTPM` option of `tximport` which adjusts the counts for differences in average transcript length between samples, so one does not need to pass a length offset matrix to the DE testing. The gene lengths are included as a column though (the median observed gene length across all samples) for normalizations such as RPKM which can be useful for visualization.

The output folders for the individual processes will be in a directory `rnaseq_preprocess_results/`.

## Resources

The indexing is hardcoded with 30GB of RAM and 6 CPUs. That memory is required as we use the genome as mapping decoy. The quantification is set with 25GB and 6 CPUs per sample. These process resources are defined in `nextflow.config`. If schedulers such as SLURM shall be used then this should be defined in `configs/schedulers.config`. We did put a SLURM default there that allows to use `-profile slurm` which submits a job to queue called "normal" with up to 8h walltime. See the Nextflow docs for [schedulers/executors](https://www.nextflow.io/docs/latest/executor.html). 

## Software

The pipeline is fully containerized. Using either `-profile docker/singularity` will use the respective containers. Using `-profile conda` is possible but containers are recommended to ensure full software-wise consistency between runs and platforms.

## Command lines

The typical command line would be:

```
#/ here using singularity and slurm with premade index:
NXF_VER=21.10.6 \
    nextflow run main.nf \
    --publishmode 'copy' \
    --idx path/to/idx --tx2gene path/to/tx2gene.txt \
    --samplesheet path/to/samplesheet.csv \
    -with-report indexing_report.html \
    -with-trace indexing_report.log \
    -profile slurm,singularity \
    -bg > report.log

#/ including the index building:
NXF_VER=21.10.6 \
    nextflow run main.nf \
    --publishmode 'copy' \
    --genome path/to/genome.fa.gz --txtome path/to/txtome.fa.gz --gtf path/to/foo.gtf.gz \
    --samplesheet path/to/samplesheet.csv \
    -with-report report.html \
    -with-trace report.log \
    -profile slurm,singularity \
    -bg > report.log
    
```    
