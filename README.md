# rnaseq_proprocess

A Nextflow pipeline for preprocessing of RNA-seq data


## Workflow

- create a genome-decoyed index with `salmon`. This requires a reference transcriptome and genome, which can be provided as files from disk or by providing links to a server, e.g. an FTP address.
Alternatively, the user can directly provide a premade index.

- Optional trimming of the fastq files with `cutadapt`. This is not activated by default, use `--trimming` to turn it on. Default adapter sequence is the TruSeq adapter.

- Quantification with `salmon` against the created or provided index.

## Usage

**Summary of important options:**  
  
-  `--fastq`: path to fastq file (pairs), e.g. `--fastq path/to/*.fastq.gz` for single-end or `--fastq path/to/*_{1,2}.fastq.gz` for paired-end data.
The data are expected to be gzip-compressed with suffix `fastq.gz`.

- `--mode`: either "single" or "paired" depending on the sequencing type.

- `--idx`: path to a premade `salmon` index. We expect salmon v1.0.0 and later.

- `--ref_genome` and `--ref_txtome`: paths to fasta files for genome and transcriptome (expected to end `.fa.gz`) or download links, e.g. FTP address to pull it from remote.

- `--trimming`: if set then perform adapter trimming with `cutadapt`

- `--quant_libtype`: the [library type flag for salmon](https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-alignment-based-mode), by default 'A' for automatic detection.
Standard stranded libraries, e.g. NEBnext UltraII Directional Prep Kit would be 'ISR' and unstranded libraries are usually 'IU'.

All other options are intuitively named in the `nextflow.config` file.

For testing there are two profiles, `-profile test_single` and `-profile test_paired` for paired- and single-end data. [Minimal example data](https://github.com/nextflow-io/rnaseq-nf/tree/master/data/ggal) are provided in `./test`,
which are also used by the GitHub Actions testing.

There are three profiles available to automatically take care of software dependencies which are `-profile conda/docker/singularity`.  
Once Nextflow makes a stable release that supports `mamba` this profile will be added as well. The `environment.yml` contains the required software for `conda`.  
A Docker image is available for the `docker,singularity` profiles at the [Docker Hub](https://hub.docker.com/r/atpoint/rnaseq_preprocess).

For submission via SLURM on a cluster one can use `-profile slurm`. The options `--queue` and `--clusteroptions` can be used to specify a queue and any other options the scheduler may use.

## Citations

-  [nf-core project from which much code inspiration was taken from](https://nf-co.re/)

-  [Ewels et al (2020) Nature Biotechnology volume 38, pages 276–278](https://www.nature.com/articles/s41587-020-0439-x)

-  [Nextflow Docs](https://www.nextflow.io/docs/latest/index.html#)

-  [Seqera Training](https://seqera.io/training/)

-  [https://github.com/nextflow-io/rnaseq-nf](https://github.com/nextflow-io/rnaseq-nf)

-  [Merkel, D (2014). Docker: lightweight linux containers for consistent development and deployment. Linux Journal](https://dl.acm.org/doi/10.5555/2600239.2600241)

-  [Kurtzer et al (2017) Singularity: Scientific containers for mobility of compute. PLoS ONE](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459)

-  [Patro et al (2017) Salmon: fast and bias-aware quantification of transcript expression using dual-phase inference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5600148/)
