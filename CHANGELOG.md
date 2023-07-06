# Changelog

## v2.5.2
- check in `tximport` process whether there is a mismatch between tx2gene file and quant.sf identifiers
that can be solved by using either of the `ignoreTxVersion` or `ignoreAfterBar` arguments of `tximport()`,
and if so, do it automatically. This is relevant for Ensembl-based annotations where the GTF file often does not contain
transcript versions (like ENSTXXXX.[digit]) but the transcriptome fasta headers do. This then creates a mismatch between
tx2gene and the quant.sf transcript names which is annoying. This check will compensate for it.

## v2.5.1
- update container

## v2.5
- adapted the idea of the meta-map from nf-core for module inputs
- merge samples with multiple fastq files prior to fastqc/quantification
- robustified reading from samplesheet
- modules that use fastq now read a `[meta, reads]` map that includes the information on
libtype and single-end/paired-end 
- added a trimmer/seqtk module
- all relevant modules now return a tuple with software versions and command lines that
a dedicated process collects and publishes in a summary file at the end of the pipeline
- updated container to v1.6.1
- added CITATIONS.md
- check if files in samplesheet exist

## v2.4
- moved samplesheet validation to a process
- updated container
- report command lines and software versions for salmon and tximport, and emit to publishDir

## v2.3
- use `stageAs` in fastqc and quant modules to deal with situation where several fastq files
per sample have the same name (name collision)

## v2.2a
- add `--skip_tximport` flag

## v2.2
- updated software versions, see `CONTAINERLOG.md`
- added a `schema.version` option to display the current pipeline version

## v2.1
- add native samplesheet validation, checking for existance of the fastq files, that technical replicates have same libtype, no mix of single/paired-end data and that header line of samplesheet is correct

## v2.0
- updated with [nf_blank](https://github.com/ATpoint/nf_blank) template
  - params validation with params listed in `schema.nf`
  - intro message with summary of all params upon startup of pipeline
  - validate minimal NF version
- read fastq files via samplesheet
- detect and automerge technical/lane replicates
- added fastqc and multiqc
- added a CI test for singularity
- add tximport

## phd_project_version
- initial release
