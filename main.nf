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
println "Documentation: https://github.com/ATpoint/rnaseq_preprocess/README.md"
println "Author:        Alexander Toenges (@ATpoint)"
println "Runname:       $workflow.runName"
println "Profile:       $workflow.profile"
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

    if(!(new File(params.samplesheet)).exists()){
        println "\u001B[31m$longline"
        println "[VALIDATION ERROR]"
        println "The samplesheet does not exist!"
        println "$longline\u001B[0m"
        System.exit(1)
    } else {
        new File(params.samplesheet).withReader { line = it.readLine() } 
        if(line != 'sample,r1,r2,libtype'){
            println "\u001B[31m$longline"
            println "[VALIDATION ERROR]"
            println "The samplesheet header must be:"
            println "sample,r1,r2,libtype -- so comma-separated, no whitespaces and exactly these titles!"
            println "$longline\u001B[0m"
            System.exit(1)
        }
    }

    // Read samplesheet and replace relative paths by absolute ones
    ch_samplesheet_initial = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header:true)
        .map {

            sample = it['sample']
            
            r1 = it['r1']
                    .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                    .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                    .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

            r2 = it['r2']
                .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

            lt = it['libtype']

            tuple(sample, r1, r2, lt)      

        }
        
        // validate:
        not_existy = [:]
        ch_samplesheet_initial.subscribe onNext: { 
            
            // validate r1 exists
            if(!(new File(it[1]).exists())) not_existy[it[0]] = it[1]

            // validate r2 exists:
            if(r2!=''){
                if(!(new File(it[2]).exists())) not_existy[it[0]] = it[2]
            }
            

        }, onComplete: { 

            // throw error and exit if any of the fastq files does not exist -- list which
            if(not_existy.size() > 0){
                println "\u001B[31m$longline"
                println "[VALIDATION ERROR]"
                println "These fastq paths do not exist as files:"
                println not_existy.toString().replaceAll("\\[|\\]", "").replaceAll(", ", "\n").replaceAll(":", " => ")
                println "$longline\u001B[0m"
                System.exit(1)
            } 

        }

    ch_samplesheet = ch_samplesheet_initial.groupTuple(by:0)

    // validate that technical reps do have the same libtype
    more_one_libtype = []
    ch_samplesheet.subscribe onNext: { 
        itl = it[3].unique()
        if(itl.size() > 1){
            more_one_libtype.add(it[0])
        }
    }, onComplete: { 
        if(more_one_libtype.size() > 0){
            println "\u001B[31m$longline"
            println "[VALIDATION ERROR]"
            println "These sample have > 1 fastq file (pair) but different libtypes:"
            println more_one_libtype.toString().replaceAll("\\[|\\]", "").replaceAll(", ", "\n").replaceAll(":", " => ")
            println "$longline\u001B[0m"
            System.exit(1)
        }
    }     

    // validate that there is no mix of single-end paired-end libraries
    end_mixes = []
    ch_samplesheet.subscribe onNext: { 
        
        rr1 = it[1] - ''      
        rr2 = it[2] - ''

        if(rr2.size() > 0){

            if(rr1.size() != rr2.size()) end_mixes.add(it[0])

        }
        

    }, onComplete: { 
        if(end_mixes.size() > 0){
            println "\u001B[31m$longline"
            println "[VALIDATION ERROR]"
            println "These sample have a mix of single-and paired-end libraries which is not supported!"
            println end_mixes.toString().replaceAll("\\[|\\]", "").replaceAll(", ", "\n").replaceAll(":", " => ")
            println "$longline\u001B[0m"
            System.exit(1)
        }
    }  


} // end of samplesheet validation   

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

def skip_quant = params.only_idx==true ? true : false
def skip_multiqc = params.only_idx==true ? true : false
def run_idx = ConvertBool2String(params.idx)=='' ? true : false
def only_fqc = params.only_fastqc==true ? true : false
def skip_tximport = params.skip_tximport==true ? true : false

workflow {

    //------------------------------------------
    // Indexing
    //------------------------------------------
    if(params.only_idx==true) {
        
        IDX()
        use_idx     = IDX.out.idx
        use_tx2gene = IDX.out.tx2gene

    } else {

        if(run_idx){

            if(run_idx==true) {

                if(only_fqc==false){
                
                    IDX()
                    use_idx     = IDX.out.idx
                    use_tx2gene = IDX.out.tx2gene

                }

            }

        } else {

            use_idx     = params.idx
            use_tx2gene = params.tx2gene

        }

    }

    //------------------------------------------
    // Fastqc
    //------------------------------------------
    if(params.skip_fastqc==false && params.only_idx==false) FASTQC(ch_samplesheet)
    
    //------------------------------------------
    // Quant
    //------------------------------------------
    if(params.only_fastqc==false && skip_quant==false){

        QUANT(ch_samplesheet, use_idx, use_tx2gene)

        quant_only_quant = QUANT.out.quant.map { k -> k[1] }

        if(!skip_tximport) { TXIMPORT(quant_only_quant.collect(), params.tximport_name, use_tx2gene) }

    }

    //------------------------------------------
    // Multiqc
    //------------------------------------------
    
    if(skip_multiqc==false){
    
        if(params.skip_fastqc == false && skip_quant == false && params.only_fastqc == false) combined_channel = FASTQC.out.fastqc.concat(quant_only_quant)
        if((params.skip_fastqc == false && skip_quant == true) || (params.only_fastqc == true))  combined_channel = FASTQC.out.fastqc
        if(params.skip_fastqc == true  && skip_quant == false) combined_channel = quant_only_quant        

        if(skip_quant == false || params.skip_fastqc == false) MULTIQC(combined_channel.collect())

    }

    def od = params.outdir
    workflow.onComplete {
        Date date2 = new Date()
        String datePart2 = date2.format("yyyy-dd-MM -- ")
        String timePart2 = date2.format("HH:mm:ss")
        def end_date = datePart2 + timePart2
        println ""
        println "\u001B[33m========================================================================================================================="
        println "Pipeline completed!"
        println "End: $end_date"
        println "Results are in:"
        println od
        println "=========================================================================================================================\u001B[0m"
        println ""
    }

}