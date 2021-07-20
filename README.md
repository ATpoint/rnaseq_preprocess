# rnaseq_proprocess

A Nextflow pipeline for preprocessing of RNA-seq data


## Workflow

- create a genome-decoyed index with `salmon`. This requires a reference transcriptome and genome, which can be provided as files from disk or by providing links to a server, e.g. an FTP address.
Alternatively, the user can directly provide a premade index.

- trimming of the fastq files with `cutadapt`. This is not activated by default, use `--trim` to turn it on. Default adapter sequence is the TruSeq adapter.

- quantification with `salmon` against the created or provided index.

- aggregation of the transcript abundance estimates to the gene level with `tximport`

All individual steps can optionally be turned off.

## Usage

A typical workflow would be to first create an index (and store it somewhere) and then quantify the data. We use the provided test data here with minimal resources.
Threads and memory should be increased for full-size data. One should probably allow up to 40GB of RAM for indexing a genome-decoyed mouse transcriptome.

First, create an index and `move` (rather than symlink which is the Nextflow default) it into the output directory:

```bash

#/ pull the repository:
git clone https://github.com/ATpoint/rnaseq_preprocess && cd rnaseq_preprocess

nextflow run main.nf \
    --ref_txtome $(realpath ./test/txtome.fa.gz) \
    --ref_genome $(realpath ./test/genome.fa.gz) \
    --ref_gtf $(realpath ./test/annot.gtf.gz) \
    --idx_threads 1 \
    --idx_mem '1.GB' \
    --outdir $(realpath ./test/) \
    --idx_name 'salmonIndex' \
    --publishdir_mode 'move' \
    --skip_quant    

```

The `test` folder now contains a subfolder `index` which stores:
- the `gentrome.fa.gz` which is the merge of transcriptome and full (decoy) genome. 
- a `tx2gene.txt` map which connects transcripts to its genes, can be used for the `tximport` step
- the index itself, here a folder called `salmonIndex`

Now quantify some paired-end fastq files, trimming the TruSeq adapter. 

```bash

nextflow run main.nf \
    --idx $(realpath ./test/index/salmonIndex) \
    --mode 'paired' \
    --trim \
    --adapter_seq 'AGATCGGAAGAGC' \
    --trim_threads 1 \
    --quant_threads 1 \
    --quant_mem '1.GB' \
    --quant_libtype 'ISR' \
    --quant_additional '\--numGibbsSamples 64' \
    --tx2gene $(realpath ./test/index/tx2gene.txt) \
    --fastq "$(realpath ./test)"'/*_{1,2}.fastq.gz' \
    --outdir $(realpath ./test)

```

The single/double quotes in `--fastq "$(realpath ./test)"'/*_{1,2}.fastq.gz'` are intended as the first part is for bash to find the full path to the data directory,
while the second part (in single quotes) is evaluated by Nextflow itself and should not be expanded by bash.

The `test` (`--outdir`) folder now contains:
- `fastq_trimmed` with the trimmed fastq files which can be deleted. There is unfortunately by now no in-built function in Nextflow to declare temporary files and folders such as this one for automated deletion.
- `salmon_quants` which stores the salmon quantification directories
- `tximport` which stores the gene level counts and average transcript lengths per sample (and inferential replicates if that option was used in `salmon`)

**Summary of important options:**  
  
-  `--fastq`: path to fastq file (pairs), e.g. `--fastq path/to/*.fastq.gz` for single-end or `--fastq path/to/*_{1,2}.fastq.gz` for paired-end data.
The data are expected to be gzip-compressed with suffix `fastq.gz`.

- `--mode`: either "single" or "paired" depending on the sequencing type.

- `--idx`: path to a premade `salmon` index. We expect salmon v1.0.0 and later.

- `--*_mem` arguments for the idx, trim, quant and tximport processes. The defaults should do for trim. For indexing a full mouse or human genome one should probably allow 40GB (for human maybe more).
For quantification 20GB should do. The tximport default should be fine unless one produced many inferential replicates. This is not tested though.

- `--ref_genome`, `--ref_txtome` and `--ref_gtf`: paths to files for genome, transcriptome (expected to end `.fa.gz`) and GTF (`gtf.gz`) or download links, e.g. FTP address to pull it from remote.
Pulling from remote in our experience is a bit unstable, therefore providing links to on-disk data is preferred.

- `--trim`: if set then perform adapter trimming with `cutadapt` using either the default (TruSeq) sequence or the one provided with `--adapter_seq`.

- `--quant_libtype`: the [library type flag for salmon](https://salmon.readthedocs.io/en/latest/salmon.html#quantifying-in-alignment-based-mode), by default 'A' for automatic detection.
Standard stranded libraries, e.g. NEBnext UltraII Directional Prep Kit would be 'ISR' and unstranded libraries are usually 'IU'.

All other options are intuitively named in the `nextflow.config` file. Some processes have an `additional` parameter to pass further arguments to the respective tool. If so the first parameter must be escaped to be parsed properly, e.g.:  
`--quant_additional '\--noLengthCorrection --numGibbsSamples 64'`

For testing there are two profiles, `-profile test_single` and `-profile test_paired` for paired- and single-end data. [Minimal example data](https://github.com/nextflow-io/rnaseq-nf/tree/master/data/ggal) are provided in `./test`,
which are also used by the GitHub Actions testing.

## Software

There are three profiles available to automatically take care of software dependencies which are `-profile conda/docker/singularity`.  
The `environment.yml` contains the required software for `conda`. Nextflow will create an enviroment based on this if this profile is set. 
A Docker image is available for the `docker,singularity` profiles at the [Docker Hub](https://hub.docker.com/r/atpoint/rnaseq_preprocess).
If none of the above profiles is specified then the software is expected locally in `PATH`. One could also use `mamba` to manually install the software, e.g.
via `mamba env update --file environment.yml`. By now there is no stable nextflow release that supports mamba out of the box. Once this changes a mamba profile will be added.

For submission via SLURM on a cluster one can use `-profile slurm`. The options `--queue` and `--clusteroptions` can be used to specify a queue and any other options the scheduler may use.

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
