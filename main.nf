#! /usr/bin/env nextflow

nextflow.enable.dsl=2

//------------------------------------------------------------------------
// Intro message
//------------------------------------------------------------------------

def longline="========================================================================================================================="

Date date = new Date()
String datePart = date.format("yyyy-dd-MM -- ")
String timePart = date.format("HH:mm:ss")
def start_date = datePart + timePart

println ""
println "\u001B[33m$longline"
println "Pipeline:      rnaseq_preprocess"
println "GitHub:        https://github.com/ATpoint/rnaseq_preprocess"
println "Documentation: https://github.com/ATpoint/rnaseq_preprocess/README.md"
println "Author:        Alexander Toenges (@ATpoint)"
println "Runname:       $workflow.runName"
println "Profile:       $workflow.profile"
println "Start:         $start_date"
println "$longline\u001B[0m"

//------------------------------------------------------------------------
// Validate input params via schema.nf
//------------------------------------------------------------------------

 evaluate(new File("${baseDir}/functions/validate_schema_params.nf"))

//------------------------------------------------------------------------
// Load the modules and pass params
//------------------------------------------------------------------------

include{ Idx }          from './modules/index'      addParams(  outdir:         params.idx_dir,
                                                                publishmode:    params.publishmode,
                                                                additional:     params.idx_additional)

include { Tx2Gene }     from './modules/tx2gene'    addParams(  outdir:         params.idx_dir,
                                                                publishmode:    params.publishmode)
                                                             
include{ FastQC }       from './modules/fastqc'     addParams(  outdir:         params.fastqc_dir,
                                                                publishmode:    params.publishmode,
                                                                additional:     params.fastqc_additional)

include{ Quant }        from './modules/quant'      addParams(  outdir:         params.quant_dir,
                                                                publishmode:    params.publishmode,
                                                                additional:     params.quant_additional)

include{ Tximport }     from './modules/tximport'   addParams(  outdir:         params.tximport_dir,
                                                                publishmode:    params.publishmode)

include{ MultiQC }      from './modules/multiqc'     addParams( outdir:         params.multiqc_dir,
                                                                publishmode:    params.publishmode,
                                                                additional:     params.multiqc_additional)

//------------------------------------------------------------------------
// Validate samplesheet
//------------------------------------------------------------------------

if(!params.only_idx){

    // Validate that fastq files in samplesheet exist as files on disk

    fastq_no_exist  = [:]
    libtype_error   = [:]
    
    if(!(new File(params.samplesheet)).exists()){
        println "\u001B[31m$longline"
        println "[VALIDATION ERROR]"
        println "The samplesheet does not exist!"
        println "$longline\u001B[0m"
        System.exit(1)
    }

    // Read samplesheet and replace relative paths by absolute ones
    ch_samplesheet = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map{ r -> 

            sample = r['sample']

            r1 = r['r1']
                    .toString()
                    .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                    .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                    .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))    

            r2 = r['r2']
                    .toString()
                    .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                    .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                    .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

            lt = r['libtype']                                                         

            tuple(sample, r1, r2, lt)                    

    }.groupTuple(by: 0)
     .map { rr -> 
                libtype = rr[3]
                if(libtype.unique().size() == 1) lt = libtype.unique()
                if(libtype.unique().size() > 1) lt = "error"

            tuple(rr[0], rr[1], rr[2], lt)

     }

}    

//------------------------------------------------------------------------      
// Define subworkflows
//------------------------------------------------------------------------

def ConvertBool2String(indata='') {
    
    if(indata instanceof Boolean){
        return ''
    } else {
        return indata
    }
}

workflow IDX {

    main:
        Idx(params.txtome, params.genome, params.idx_name)

        Tx2Gene(params.gtf)

        this_idx     = Idx.out.idx 
        this_tx2gene = Tx2Gene.out.tx2gene
            
    emit:
        idx     = this_idx
        tx2gene = this_tx2gene    

}

workflow FASTQC {

    take:
        samplesheet

    main:

        new_samplesheet = samplesheet.map { k ->

            // maybe one day there will be optional inputs for DSL2 ...
            r2 = k[2].toString() == "['']" ? "/" : k[2]
            tuple(k[0], k[1], r2)

        }

        FastQC(new_samplesheet)

    emit:
        fastqc = FastQC.out.zip


}

workflow QUANT {

    take:
        samplesheet
        idx
        tx2gene
        
    main:

       new_samplesheet = samplesheet.map { k ->

            // maybe one day there will be optional inputs for DSL2 ...
            r2 = k[2].toString() == "['']" ? "/" : k[2]
            tuple(k[0], k[1], r2, k[3])

        }

        Quant(new_samplesheet, idx, tx2gene)

    emit:
        quant   = Quant.out.quants
        tx2gene = Quant.out.tx2gene

}

workflow TXIMPORT {

    take:
        salmons
        outname
        tx2gene

    main:
        Tximport(salmons, outname, tx2gene)
}

workflow MULTIQC {

    take:
        paths
        
    main:

       MultiQC(paths)

}

//------------------------------------------------------------------------      
// Define and run the main workflow
//------------------------------------------------------------------------

if(params.only_idx) {
    skip_fastqc = true
    skip_quant = true
    skip_multiqc = true
} else {
    skip_fastqc = false
    skip_quant = false
    skip_multiqc = false
}

if(params.only_fastqc==true){
    only_fastqc = true
} else {
    only_fastqc = false
}

if(params.skip_fastqc==true){
    skip_fastqc == true
} else {
    skip_fastqc == false
}
        
workflow {

    //------------------------------------------
    // Indexing
    //------------------------------------------
    if(params.only_idx) {
        
        IDX()
        use_idx     = IDX.out.idx
        use_tx2gene = IDX.out.tx2gene

    } else {

        if(ConvertBool2String(params.idx)==''){

            if(only_fastqc==false) {
                IDX()
                use_idx     = IDX.out.idx
                use_tx2gene = IDX.out.tx2gene
            }

        } else {

            use_idx     = params.idx
            use_tx2gene = params.tx2gene

        }

    }

    //------------------------------------------
    // Fastqc
    //------------------------------------------
    if(skip_fastqc==false) FASTQC(ch_samplesheet)
    
    //------------------------------------------
    // Quant
    //------------------------------------------
    if(only_fastqc == false && skip_quant == false){

        QUANT(ch_samplesheet, use_idx, use_tx2gene)

        quant_only_quant = QUANT.out.quant.map { k -> k[1] }

        TXIMPORT(quant_only_quant.collect(), params.tximport_name, use_tx2gene)

    }

    //------------------------------------------
    // Multiqc
    //------------------------------------------
    if(skip_multiqc==false){
    
        if(skip_fastqc == false && skip_quant == false && only_fastqc == false) combined_channel = FASTQC.out.fastqc.concat(quant_only_quant)
        if((skip_fastqc == false && skip_quant == true) || (only_fastqc == true))  combined_channel = FASTQC.out.fastqc
        if(skip_fastqc == true  && skip_quant == false) combined_channel = quant_only_quant        

        if(skip_quant == false || skip_fastqc == false) MULTIQC(combined_channel.collect())

    }

}

