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
println ""
println "Pipeline:      rnaseq_preprocess"
println "GitHub:        https://github.com/ATpoint/rnaseq_preprocess"
println "Revision:      $workflow.revision"
println "Commit:        $workflow.commitId"
println "Documentation: https://github.com/ATpoint/rnaseq_preprocess/README.md"
println "Author:        Alexander Toenges (@ATpoint)"
println "Runname:       $workflow.runName"
println "Profile:       $workflow.profile"
println "Command line:  $workflow.commandLine"
println "Start:         $start_date"
println ""
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

include{ ValidateSamplesheet }  from './modules/validatesamplesheet'

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

include{ CommandLines } from './modules/commandline' addParams( outdir:         params.pipedir,
                                                                publishmode:    params.publishmode)
                                                              

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
        cl      = Idx.out.commandlines
        vs      = Idx.out.versions

}

workflow VALIDATESSAMPLESHEET {

    take: 
        samplesheet_in

    main:
        ValidateSamplesheet(samplesheet_in)

    sout = ValidateSamplesheet.out.samplesheet
           .splitCsv(header:true)
           .map {

                sample = it['sample']
                
                r1 = it['r1']
                        .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                        .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                        .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

                rx = it['r2']
                    .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                    .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                    .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

                r2 = rx.toString()=='' ? '/' : rx

                lt = it['libtype']

                tuple(sample, r1, r2, lt)      
                
            }
            .groupTuple(by:0)

    emit:
        samplesheet = sout

}

workflow FASTQC {

    take:
        samplesheet

    main:

        FastQC(samplesheet)

    emit:
        fastqc = FastQC.out.zip


}

workflow QUANT {

    take:
        samplesheet
        idx
        tx2gene
        
    main:

        Quant(samplesheet, idx, tx2gene)

    emit:
        quant   = Quant.out.quants
        tx2gene = Quant.out.tx2gene
        cl      = Quant.out.commandlines
        vs      = Quant.out.versions

}

workflow TXIMPORT {

    take:
        salmons
        tx2gene

    main:
        Tximport(salmons, tx2gene)

    emit:
        vs      = Tximport.out.versions

}

workflow MULTIQC {

    take:
        paths
        
    main:

       MultiQC(paths)

}

workflow EVERYTHING {

    if(params.only_idx==true) {
        
        IDX()
        use_idx     = IDX.out.idx
        use_tx2gene = IDX.out.tx2gene
        CommandLines(IDX.out.cl.collect(), IDX.out.vs.collect())

    } else {

        VALIDATESSAMPLESHEET(params.samplesheet)

        if(params.skip_fastqc!=true){
            FASTQC(VALIDATESSAMPLESHEET.out.samplesheet)
            fastqc_for_multiqc = FASTQC.out
        } else fastqc_for_multiqc = Channel.empty()

        if(params.only_fastqc==false & params.only_idx==false){

            use_idx     = params.idx
            use_tx2gene = params.tx2gene

            QUANT(VALIDATESSAMPLESHEET.out.samplesheet, use_idx, use_tx2gene)

            quant_for_multiqc = QUANT.out.quant

            if(!params.skip_tximport) { 
                
                TXIMPORT(QUANT.out.quant.collect(), use_tx2gene)

                CommandLines(QUANT.out.cl.collect(),
                             QUANT.out.vs.concat(TXIMPORT.out.vs).collect())


            } else {

                CommandLines(QUANT.out.cl.collect(), QUANT.out.vs.collect())

            }

        } else quant_for_multiqc = Channel.empty()

        MULTIQC(fastqc_for_multiqc.concat(quant_for_multiqc).collect())

    }

}

workflow { EVERYTHING() }

def od = params.outdir
workflow.onComplete {
    Date date2 = new Date()
    String datePart2 = date2.format("yyyy-dd-MM -- ")
    String timePart2 = date2.format("HH:mm:ss")
    def end_date = datePart2 + timePart2
    println ""
    println "\u001B[33m$longline"
    println "Pipeline completed!"
    println "End: $end_date"
    println "Results are in:"
    println od
    println "$longline\u001B[0m"
    println ""
}