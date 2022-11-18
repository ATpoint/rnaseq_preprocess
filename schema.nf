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
schema.version        = [value: 'v2.3', type: 'string']
schema.min_nf_version = [value: '21.10.6', type: 'string', mandatory: true, allowed: '']
schema.publishmode    = [value: 'copy', type: 'string', mandatory: true, allowed:['symlink', 'rellink', 'link', 'copy', 'copyNoFollow', 'move']]
overall_outdir        = "$launchDir/rnaseq_preprocess_results/"
schema.outdir         = [value: overall_outdir, type: 'string', mandatory: true]
schema.pipedir        = [value: "${overall_outdir}/pipeline/", type: 'string']


// indexing:
schema.title2          = [title: 'INDEXING OPTIONS']
schema.idx             = [value: '', type: 'string']
schema.tx2gene         = [value: '', type: 'string']
schema.txtome          = [value: '', type: 'string', pattern: /.*\.gz$/]
schema.genome          = [value: '', type: 'string', pattern: /.*\.gz$/]
schema.gtf             = [value: '', type: 'string', pattern: /.*\.gz$/]
schema.transcript_id   = [value: 'transcript_id', type: 'string']
schema.transcript_name = [value: 'transcript_name', type: 'string']
schema.gene_id         = [value: 'gene_id', type: 'string']
schema.gene_name       = [value: 'gene_name', type: 'string']
schema.gene_type       = [value: 'gene_type', type: 'string']
schema.idx_name        = [value: "idx", type: 'string']
schema.idx_dir         = [value: "${overall_outdir}/salmonIndex/", type: 'string']
schema.idx_additional  = [value: '--gencode', type: 'string']
schema.only_idx        = [value: false, type: 'logical']

// fastqc/multiqc
schema.title3             = [title: 'FASTQC/MULTIQC OPTIONS']
schema.skip_fastqc        = [value: false, type: 'logical']
schema.fastqc_dir         = [value: "${overall_outdir}/fastqc/", type: 'string']
schema.multiqc_dir        = [value: "${overall_outdir}/multiqc/", type: 'string']
schema.fastqc_additional  = [value: '', type: 'string']
schema.multiqc_additional = [value: '', type: 'string']
schema.only_fastqc        = [value: false, type: 'logical']

// samplesheet and quantification
schema.title3             = [title: 'QUANTIFICATION OPTIONS']
schema.samplesheet        = [value: '', type: 'string', pattern: /.*\.csv$/]
schema.quant_dir          = [value: "${overall_outdir}/salmonQuant/", type: 'string']
schema.quant_additional   = [value: '--gcBias --seqBias --posBias', type: 'string']
schema.tximport_dir       = [value: "${overall_outdir}/tximport/", type: 'string']
schema.skip_tximport      = [value: false, type: 'logical']

// related to the container/environment for the R/Bioconductor part of this workflow
schema.title4         = [title: 'CONTAINER/CONDA OPTIONS']
schema.container      = [value:'atpoint/rnaseq_preprocess:v1.5.0', type:'string', mandatory:true]
schema.environment    = [value: "$baseDir/environment.yml", type:'string', mandatory: true ]

// --------------------------------------------------------------------------------------------------------------

return schema // don't change this line
