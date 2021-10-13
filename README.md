# rnaseq_proprocess

A Nextflow pipeline for preprocessing of RNA-seq data.

<br>

![CI](https://github.com/ATpoint/rnaseq_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Docker builds](https://img.shields.io/docker/automated/atpoint/rnaseq_preprocess]


<br>

The pipelines quantifies RNA-seq data in fastq format (fastq.gz), optionally after trimming, against a reference using [salmon](https://github.com/COMBINE-lab/salmon). It then summarizes the transcript abundance estimates to the gene level with [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html). 

The user can run the pipeline with software being installed locally (default), or via either conda, Docker and Singularity by using `-profile conda/docker/singularity`. This will then either create the conda environment based on the `environment.yml` in this repository or pull [an image with the required software](https://hub.docker.com/r/atpoint/rnaseq_preprocess/tags?page=1&ordering=last_updated) for Docker/Singularity automatically.

Defaults for the parameters can be found in `nextflow.config`. 

- params.publishdir_mode : publish mode for results, see [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir)
- params.mode            : paired/single for paired-end and single-end data. File names must be `*.fastq.gz` for single-end and `*_1/2.fastq.gz` for pared-end data.
- params.ref_txtome      : path to a reference transcriptome in (fa.gz) to build index from
- params.ref_genome      : path to reference genome (fa.gz) to use as decoy during indexing
- params.ref_gtf         : path to reference GTF to parse tx2gene information from (linking transcripts to its parent gene
- params.fastq           : path to input data (fastq.gz)
- params.idx             : path to the index folder in case a premade index shall be used, this then skips the indexing step
- params.idx_name        : name of the index being created
- params.idx_threads     : threads for indexing
- params.idx_mem         : memory for indexing
- params.idx_additional  : additional command line parameters for `salmon index`
- params.idx_dir         : outdir for index (`params.idx.name` is the index folder generated in this directory)
- params.transcript_id   : name of the GTF column that stores transcript ID
- params.transcript_name : name of the GTF column that stores transcript name
- params.gene_id         : name of the GTF column that stores gene ID
- params.gene_name       : name of the GTF column that stores gene name
- params.gene_type       : name of the GTF column that stores gene type
- params.trim            : logical, whether to trim fastq files for adapters
- params.trim_additional : additional command line parameters for `cutadapt`
- params.trim_dir        : output directory for trimmed fastq files
- params.adapter_seq     : adapter sequence to trim (currently only a single sequence is supported, typically this would be TruSeq universal sequence)
- params.trim_threads    : threads for trimming
- params.trim_mem        : memory for trimming
- params.skip_quant      : logical, whether to skip quantification step
- params.quant_threads   : threads for quantification
- params.quant_mem       : memory for qualtification
- params.quant_additional: additional command line arguments for `salmon quant`
- params.quant_libtype   : [library type](https://salmon.readthedocs.io/en/latest/library_type.html) for quantification
- params.quant_dir       : output directory for quantification results
- params.skip_tximport   : logical, whether to skip tximport step
- params.tx2gene         : path to tx2gene file (created during indexing)
- params.tximport_dir    : output directory for tximport results
- params.tximport_mem    : memory for tximport
- params.queue           : SLURM queue name to submit to if using SLURM via `-profile slurm`
- params.clusteroptions  : additional options for cluster submission

<br>
Example launch command for a standard paired-end dataset, with index already created, using Singularity, and submitting via SLURM on a 'normal' queue,
using up to 16 threads per quantification job:
<br>

```bash
module load Singularity && \
NXF_HOME=$(realpath ./NXF_HOME/) NXF_VER=21.04.1 \
    nextflow run atpoint/rnaseq_preprocess -r 8ba0a6bf7866d2272c44ae33613373134af41ae4 \
        -with-trace \
        -profile singularity,slurm \
        --queue 'normal' \
        --fastq "$(pwd)"'/*{1,2}.fastq.gz' \
        --idx '/path/to/idx_folder/' \
        --trim --trim_threads 8 \
        --quant_threads 16 --quant_libtype 'A' --quant_additional '\--gcBias --seqBias' \
        --tx2gene '/path/to/tx2gene.txt'      
```

For testing there are small example data in the `test` folder, and three profiles to use them (`-profile test_paired/test_single/test_idx`). These are used for the GitHub Actions [CI testing](https://github.com/ATpoint/rnaseq_preprocess/blob/main/.github/workflows/CI.yml).

## Citations

-  [nf-core project from which much code inspiration was taken from](https://nf-co.re/)

-  [Ewels et al (2020) The nf-core framework for community-curated bioinformatics pipelines](https://www.nature.com/articles/s41587-020-0439-x)

-  [The Nextflow Documentation](https://www.nextflow.io/docs/latest/index.html#)

-  [Seqera Training](https://seqera.io/training/)

-  [https://github.com/nextflow-io/rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf)

-  [Merkel, D (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal](https://dl.acm.org/doi/10.5555/2600239.2600241)

-  [Kurtzer et al (2017) Singularity: Scientific containers for mobility of compute. PLoS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459)

-  [Patro et al (2017) Salmon: fast and bias-aware quantification of transcript expression using dual-phase inference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/)

-  [Soneson et al (2015) Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences](https://f1000research.com/articles/4-1521/v2)
