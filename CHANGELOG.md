# Changelog

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
