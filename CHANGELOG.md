# Changelog

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