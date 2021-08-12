# rnaseq_proprocess

A Nextflow pipeline for preprocessing of RNA-seq data.

<br>

![CI](https://github.com/ATpoint/rnaseq_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>


...documentation will follow...

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
        --idx '/scratch/tmp/a_toen03/project/reference_files/salmon_idx/gencode_vM25_k31/idx_k31' \
        --trim --trim_threads 8 \
        --quant_threads 16 --quant_libtype 'A' --quant_additional '\--gcBias --seqBias' \
        --tx2gene '/scratch/tmp/a_toen03/project/reference_files/salmon_idx/gencode_vM25_k31/tx2gene.txt'      
```

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
