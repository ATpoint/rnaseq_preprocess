#! /usr/bin/env nextflow

/*
 *
 *  SCHEMA DEFINITION FOR PARAMS VALIDATION
 *
 */

def Map schema = [:] // don't change this line

// --------------------------------------------------------------------------------------------------------------

// generic options:
schema.title1         = [title: 'GENERIC OPTIONS']
schema.min_nf_version = [value: '21.10.6', type: 'string', mandatory: true, allowed: '']
schema.publishmode    = [value: 'copy', type: 'string', mandatory: true, allowed:['symlink', 'rellink', 'link', 'copy', 'copyNoFollow', 'move']]
schema.outdir         = [value: "$launchDir/rnaseq_preprocess_results/", type: 'string', mandatory: true]
schema.pipe_dir        = [value: "${schema.outdir['value']}/pipeline_info/", type: 'string']

// combining technical replicates into a single fastq and trimming
schema.title3            = [title: 'MERGE/FASTQ OPTIONS']
schema.merge_dir         = [value: "${schema.outdir['value']}/fastq_merged/", type: 'string']
schema.merge_keep        = [value: false, type: 'logical']
schema.trim_reads        = [value: false, type: 'logical']
schema.trim_dir          = [value: "${schema.outdir['value']}/fastq_trimmed/", type: 'string']
schema.trim_length       = [value: 75, type: 'numeric']
schema.trim_keep         = [value: false, type: 'logical']


// fastqc/multiqc
schema.title4             = [title: 'FASTQC/MULTIQC OPTIONS']
schema.skip_fastqc        = [value: false, type: 'logical']
schema.skip_multiqc       = [value: false, type: 'logical']
schema.fastqc_dir         = [value: "${schema.outdir['value']}/fastqc/", type: 'string']
schema.multiqc_dir        = [value: "${schema.outdir['value']}/multiqc/", type: 'string']
schema.only_fastqc        = [value: false, type: 'logical']

// samplesheet and quantification
schema.title5             = [title: 'QUANTIFICATION/TXIMPORT OPTIONS']
schema.samplesheet        = [value: '', type: 'string', pattern: /.*\.csv$/]
schema.idx                = [value: '', type: 'string', mandatory: true]
schema.quant_dir          = [value: "${schema.outdir['value']}/salmon_quant/", type: 'string']
schema.quant_additional   = [value: '--gcBias --seqBias --posBias', type: 'string']
schema.tx2gene            = [value: '', type: 'string']
schema.tximport_dir       = [value: "${schema.outdir['value']}/tximport/", type: 'string']
schema.skip_tximport      = [value: false, type: 'logical']

// related to the container/environment for the R/Bioconductor part of this workflow
schema.title6         = [title: 'CONTAINER/CONDA OPTIONS']
schema.container      = [value: 'atpoint/rnaseq_preprocess:v1.7.2', type: 'string', mandatory: true]
schema.environment    = [value: "$baseDir/environment.yml", type:'string', mandatory: true ]

// --------------------------------------------------------------------------------------------------------------

return schema // don't change this line
