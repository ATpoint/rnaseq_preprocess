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

def latest_sha = "git rev-parse HEAD".execute()

println ''
println '|-------------------------------------------------------------------------------------------------------------'
println ''
println "[Info] This is rnaseq_preprocess
println ''
println '|-------------------------------------------------------------------------------------------------------------'
println ''

// fastq file channel:
if(!params.skip_quant) {
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
} else ch_fastq = null

// Define the final workflow:
workflow RNASEQ {

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Indexing

    include{    SalmonIndex }       from './modules/index'      addParams(  threads:    params.idx_threads,
                                                                            mem:        params.idx_mem,
                                                                            outdir:     params.idx_dir,
                                                                            pubmode:    params.publishdir_mode,
                                                                            additional: params.idx_additional)

    
    // Either make a new index from scratch or use provided one if exists:
    if(params.idx == ''){
        
        SalmonIndex(params.ref_txtome, params.ref_genome, params.idx_name, params.ref_gtf)
        use_index = SalmonIndex.out.idx
        use_tx2gene = SalmonIndex.out.tx2gene

    } else {
        
        if(! file(params.idx).exists()){

            println("[Error] ::: --idx does not exist!")
            System.exit(1)

        } else {
            
            use_index = params.idx 

            // if no indexing from scratch then make sure tx2gene exists
            if(! new File(params.tx2gene).exists() | params.tx2gene == ''){
                println("[Error] ::: --tx2gene does not exist!")
                System.exit(1)

            } else use_tx2gene = params.tx2gene
        }



    }

    //-------------------------------------------------------------------------------------------------------------------------------//
    // Trimming

    include{    Trim        }       from './modules/trim'       addParams(  threads:    params.trim_threads,
                                                                            mem:        params.trim_mem,
                                                                            outdir:     params.trim_dir,
                                                                            pubmode:    params.publishdir_mode,
                                                                            additional: params.trim_additional)

    if(params.trim && !params.skip_quant) {
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

    if(!params.skip_quant) SalmonQuant(fq, use_index)
            
    //-------------------------------------------------------------------------------------------------------------------------------//
    // tximport

    include{    Tximport }          from './modules/tximport'   addParams(  outdir:     params.tximport_dir,
                                                                            pubmode:    params.publishdir_mode,
                                                                            mem:        params.tximport_mem)

    // only if quant was run:
    run_tximport = true

    if(params.skip_quant | params.skip_tximport) run_tximport = false

    if(run_tximport) Tximport(SalmonQuant.out.quants, use_tx2gene)

}

// Run it:
workflow { RNASEQ() }
