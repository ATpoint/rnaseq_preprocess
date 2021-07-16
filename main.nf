#! /usr/bin/env nextflow

/*  
    ---------------------------------------------------------------------------------------------------------------
    
    Workflow for RNA-seq quantification with salmon
    --||    Creation of a genome-decoyed index

    --||    optional trimming of the data

    --||    quantification with salmon

    ---------------------------------------------------------------------------------------------------------------
*/ 

nextflow.enable.dsl=2

println ''
println '|------------------------------------------------------------------'
println "[Info] This is rnaseq_preprocess version ::: $params.version"
println '|------------------------------------------------------------------'
println ''

// fastq file channel:
if(params.mode == "paired"){
    ch_fastq    = Channel
                .fromFilePairs(params.fastq, checkIfExists: true)
                .ifEmpty("No fastq files found")
} else if(params.mode == "single"){
    ch_fastq    = Channel
                .fromPath(params.fastq, checkIfExists: true)
                .ifEmpty("No fastq files found")
                .map { file -> tuple(file.simpleName, file) }
}

// Define the final workflow:
workflow RNASEQ {

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Indexing

    include{    SalmonIndex }       from './modules/index'      addParams(  threads:    params.idx_threads,
                                                                            mem:        params.idx_mem,
                                                                            outdir:     params.idx_dir,
                                                                            pubmode:    params.publishdir_mode,
                                                                            additional: params.idx_additional)

    SalmonIndex(params.ref_txtome, params.ref_genome, params.idx_name, params.ref_gtf)

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Trimming
    include{    Trim        }       from './modules/trim'       addParams(  threads:    params.trim_threads,
                                                                            mem:        params.trim_mem,
                                                                            outdir:     params.trim_dir,
                                                                            pubmode:    params.publishdir_mode,
                                                                            additional: params.trim_additional)

    if(params.trim) {
        Trim(ch_fastq)
        fq = Trim.out.trimmed
    } else fq = ch_fastq

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Quantification

    // we add some defaults for paired-end data if params.quant_additional is empty:
    q_additional = params.quant_additional
    if(params.quant_additional == '' && params.mode == 'paired'){
        q_additional = '--gcBias --seqBias'
    }

    include{    SalmonQuant }       from './modules/quant'      addParams(  threads:    params.quant_threads,
                                                                            mem:        params.quant_mem,
                                                                            outdir:     params.quant_dir,
                                                                            pubmode:    params.publishdir_mode,
                                                                            libtype:    params.quant_libtype,
                                                                            additional: q_additional)

    // either use the index that was supplied by params.idx or the one made if params.ref_txtome/ref_genome was provided
    if(!params.skip_quant){
        if(params.idx == ''){
            SalmonQuant(fq, SalmonIndex.out.idx)
        } else {
            SalmonQuant(fq, params.idx)
        }
    }
    
    //-------------------------------------------------------------------------------------------------------------------------------//
    // tximport
    include{    Tximport }          from './modules/tximport'   addParams(  outdir:     params.tximport_dir,
                                                                            pubmode:    params.publishdir_mode,
                                                                            mem:        params.tximport_mem)

    if(!params.skip_tximport){
        Tximport(SalmonQuant.out.quants, SalmonIndex.out.tx2gene)
        

    }

}

// Run it:
workflow { RNASEQ() }
