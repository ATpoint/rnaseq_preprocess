# rnaseq_preprocess

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>

## Introduction

**rnaseq_preprocess** is a Nextflow pipeline for RNA-seq quantification with `salmon`. The processing steps are `fastqc` first, then quantification, aggregation to gene level with `tximport` and a small summary report with `MultiQC`. The pipeline is containerized via Docker and Singularity (container: `atpoint/rnaseq_preprocess/v1.6.0`). Outputs can be found in `rnaseq_preprocess_results/`.

# Details

**Indexing**<br>

The indexing step must be run separately first. For this we need a reference transcriptome (gzipped), a reference genome as decoy (gzipped) and a GTF annotation file (gzipped).
The flags to provide these files are `--genome`, `--txtome` and `--gtf`. Also, for creating a transcripts to gene (tx2gene) map for gene level summarization after mapping the pipeline needs `--transcript_id`, `--transcript_name`, `--gene_id`, `--gene_name` and `--gene_type` which take the colnames of the GTF that store these attributes. The defaults assume GENCODE reference annotations, and if using them these options do not need to be provided. For this we also pass the `--gencode` flag to the `salmon` indexing process via `--idx_additional` like `--idx_additional '--gencode'`. Any other salmon indexing options can be passed via this flag as well such as non-standard kmer length, e.g. `idx_additional '--gencode -k 17'` which makes sense when quantifying very short reads such as the 2x25bp data.
The output of the indexing process will be in the folder `rnaseq_preprocess_results/salmonIdx/name` where "name" is the value passed via `--idx_name` which defaults to simply "idx".
The content of the `salmonIdx/idx/` folder is what the downstream quantification process expects via `--idx` and the created tx2gene map is expected via `--tx2gene`.
We hardcoded 30GB of RAM and 6 CPUs as resources for this.

**Command for our HPC**:
```bash
NXF_VER=21.10.6 nextflow run main.nf -profile singularity,slurm --publishmode 'copy' --only_idx \
    --genome path/to/genome.fa.gz --txtome path/to/txtome.fa.gz --gtf path/to/foo.gtf.gz \
    -with-report indexing_report.html -with-trace indexing_report.trace -bg > indexing_report.log
```

**Quantification**

The quantification requires a samplesheet with four columns such as the one in the [test folder](test/samplesheet.csv):

```bash
sample,r1,r2,libtype
sample1,$baseDir/test/sample1_1.fastq.gz,$baseDir/test/sample1_2.fastq.gz,A
sample1,$baseDir/test/test2/sample1_1.fastq.gz,$baseDir/test/test2/sample1_2.fastq.gz,A
sample2,$baseDir/test/sample2_1.fastq.gz,$baseDir/test/sample2_2.fastq.gz,A
sample3,$baseDir/test/sample3.fastq.gz,,A
sample3,$baseDir/test/test2/sample3.fastq.gz,,A
```

The first row is the header which must be `sample,r1,r2,libtype`. These four fields are:
- `sample`: The sample name, this can be any user-defined name (no whitespaces here!)    
- `r1/r2`: The paths to the fastq files belonging to this sample. If this is single-end then leave `r2` empty, the pipeline will then detect this and trigger single-end mode.
- `libtype`: The library type compatible with the `--libtype` argument of `salmon`. Default is automatic detection (`A`). For Illumina stranded paired-end data that is mostly `ISR`. See the [salmon docs](https://salmon.readthedocs.io/en/latest/library_type.html) for appropriate choices of your libraries. 

Technical/lane replicates can provided here as well by using the **same** sample name (column1) for the fastq files being a technical replicate. In the above example `sample2` has two fastq file pairs that are a technical replicate. These will be merged internally prior to quantification. There is be a pre-flight samplesheet validation process that should capture corrupted samplesheets, detect potential duplicate files etc.

For quantification, custom command line arguments can be passed to `salmon quant` (which runs the quantification) by `--quant_additional`, e.g. `--quant_additional '--numGibbsSamples 8'` to produce inferential replicates. By default the ipeline runs paired-end data with the `salmon` flags `--posBias --seqBias --gcBias` and single-end data with `--posBias --seqBias`. Use `--quant_additional ' '` to turn that off. In case of **end-tagged** data such as QuantSeq one should use `--quant_additional '\--noLengthCorrection'`.

Prior to quantification a `fastqc` check is run. Options `--only_fastqc` and `--skip_fastqc` are possible.

After quantification (unless `--skip_tximport`) the transcript abundance estimates are summarized to gene level using `tximport` with the `lengthScaledTPM` option of [tximport](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#Downstream_DGE_in_Bioconductor) so the counts are ready to go into any downstream DEG framework, and counts are already corrected for average transcript length.

The quantification is hardcoded with 25GB RAM and 6 CPUs per sample which is sufficient for human and mouse genome-decoyed indices.

Exact command lines for salmon and the software versions are emitted in the `pipeline` subfolder of the output directory.

**Command for our HPC**:
```bash
NXF_VER=21.10.6 nextflow run main.nf -profile singularity,slurm --publishmode 'copy' \
    --idx path/to/idx --tx2gene path/to/tx2gene.txt --samplesheet path/to/samplesheet.csv \
    -with-report quant_report.html -with-trace quant_report.trace -bg > quant_report.log
```