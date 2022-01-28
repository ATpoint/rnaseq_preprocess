# rnaseq_preprocess

<br>

![CI](https://github.com/ATpoint/sc_preprocess/actions/workflows/CI.yml/badge.svg)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.6-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

<br>

## Introduction

**rnaseq_preprocess** is an automated preprocessing pipeline for RNA-seq data implemented using [Nextflow](https://www.nextflow.io/) which is fully containerized to take care of all required software and ensure reproducibility. It performs an initial QC with `fastqc`, then builds a genome-decoyed index against a reference transcriptome and quantifies the fastq files against it with [salmon](https://salmon.readthedocs.io/en/latest/salmon.html) from [Rob Patro's lab](https://combine-lab.github.io/). Eventually, the individual quantifications are summarized to the gene level with `tximport` returing a single matrix of raw counts. The `fastqc` and mapping statistics are summarized with `MultiQC`. The pipeline supports automated merging of technical/lane/sequencing replicates and auto-detects whether data are single- or paired-end.

# Details

tba...

# Usage

tba...
