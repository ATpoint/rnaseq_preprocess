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

include{ Idx } from './modules/index' addParams(outdir: params.idx_dir, additional: params.idx_additional)

include{ ValidateSamplesheet } from './modules/validatesamplesheet'

include{ CatFastq } from './modules/cat_fastq' addParams(outdir: params.merge_dir, keep: params.merge_keep)

include { Tx2Gene } from './modules/tx2gene' addParams(outdir: params.idx_dir)
                                                             
include{ FastQC } from './modules/fastqc' addParams(outdir: params.fastqc_dir)

include{ Trim } from './modules/trim' addParams(outdir: params.trim_dir, keep: params.trim_keep)

include{ Quant } from './modules/quant' addParams(outdir: params.quant_dir, additional: params.quant_additional)

include{ Tximport } from './modules/tximport' addParams(outdir: params.tximport_dir)

include{ MultiQC } from './modules/multiqc' addParams(outdir: params.multiqc_dir)

include{ CommandLines } from './modules/commandline' addParams( outdir: params.pipe_dir)                                                              

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
        idx      = this_idx
        tx2gene  = this_tx2gene    
        versions = Idx.out.versions

}

workflow VALIDATESSAMPLESHEET {

    take: 
        samplesheet_unvalidated

    main:
        ValidateSamplesheet(samplesheet_unvalidated)

    samplesheet = ValidateSamplesheet.out.samplesheet
           .splitCsv(header:true)
           .map {
               
                // Samplesheet allows Nextflow variables to be used, replace by absolute path
                r1 = it['r1']
                        .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                        .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                        .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

                rx = it['r2']
                        .replaceAll('\\$baseDir|\\$\\{baseDir\\}', new String("${baseDir}/"))
                        .replaceAll('\\$launchDir|\\$\\{launchDir\\}', new String("${launchDir}/"))
                        .replaceAll('\\$projectDir|\\$\\{projectDir\\}', new String("${projectDir}/"))

                r2 = rx.toString()=='' ? '.' : rx   

                se = rx.toString()=='' ? true : false

                // meta map inspired by nf-core
                meta = [id:it['sample'], lib_type:it['libtype'], single_end: se]     
                reads = [r1: r1, r2: r2]      
                counter = [1]

                tuple(meta, reads, counter)
                
            }
            .groupTuple(by:0)
            .map{ meta, grouped_reads, counter -> [ meta, grouped_reads.flatten(), counter ] }
            .map{ meta, grouped_reads, counter -> 

                    if(meta.single_end){
                        reads = grouped_reads['r1'].flatten()
                    } else {
                        reads = [grouped_reads['r1'].flatten(), grouped_reads['r2'].flatten()].flatten()
                    }

                    // counter is > 1 if a sample has more than one fastq file per read, requiring a merge
                    // before actual fastq processing
                    [meta, reads, counter.size()]

            }

    emit:
        samplesheet = samplesheet

}

workflow FASTQC {

    take:
        samplesheet

    main:

        FastQC(samplesheet)

    emit:
        fastqc = FastQC.out.zip
        versions = FastQC.out.versions


}

workflow TRIM {

    take:
        samplesheet

    main:

        Trim(samplesheet)

    emit:
        fastq_tuple = Trim.out.fastq_tuple
        versions = Trim.out.versions

}

workflow QUANT {

    take:
        samplesheet
        idx
        
    main:

        Quant(samplesheet, idx)

    emit:
        quant   = Quant.out.quant
        versions = Quant.out.versions

}

workflow TXIMPORT {

    take:
        salmons
        tx2gene

    main:
        Tximport(salmons, tx2gene)

    emit:
        versions = Tximport.out.versions

}

workflow MULTIQC {

    take:
        paths
        
    main:

       MultiQC(paths)

}

workflow RNASEQ_PREPROCESS {

    if(params.only_idx==true) {
        
        IDX()
        use_idx     = IDX.out.idx
        use_tx2gene = IDX.out.tx2gene
        idx_versions = IDX.out.versions

        cat_versions = Channel.empty()
        fastqc_versions = Channel.empty()
        trim_versions = Channel.empty()
        quant_versions = Channel.empty()
        tximport_versions = Channel.empty()

    } else {

        idx_versions = Channel.empty()
        use_idx = params.idx
        use_tx2gene = params.tx2gene

        // ----------------------------------------------------------------------------------------
        // Validate the provided samplesheet and merge fastq if necessary
        // ----------------------------------------------------------------------------------------

        VALIDATESSAMPLESHEET(params.samplesheet)

        // Samples with > 1 fastq per read
        VALIDATESSAMPLESHEET.out.samplesheet
        .map {meta, reads, counter -> 
        
            if(counter>1) [meta, reads]

        }.set { ch_needMerge }

        CatFastq(ch_needMerge)
        ch_merged = CatFastq.out.fastq_tuple
        cat_versions = CatFastq.out.versions  

        // Samples with 1 fastq per read
        VALIDATESSAMPLESHEET.out.samplesheet
        .map {meta, reads, counter -> 
            
            if(counter==1) [meta, reads]

        }.set { ch_noMerge }

        // This channel is now [meta, reads] and can go in all downstream processes that require fastq
        ch_fastq = ch_noMerge.concat(ch_merged)

        // ----------------------------------------------------------------------------------------
        // Fastqc
        // ----------------------------------------------------------------------------------------

        if(!params.skip_fastqc){

            FASTQC(ch_fastq)
            fastqc_for_multiqc = FASTQC.out.fastqc
            fastqc_versions = FASTQC.out.versions

        } else {

            fastqc_for_multiqc = Channel.empty()
            fastqc_versions = Channel.empty()

        }

        // ----------------------------------------------------------------------------------------
        // Trim
        // ----------------------------------------------------------------------------------------
        if(params.trim_reads){

            TRIM(ch_fastq)
            reads_for_quant = TRIM.out.fastq_tuple
            trim_versions = TRIM.out.versions

        } else {

            reads_for_quant = ch_fastq
            trim_versions = Chennel.empty()

        }

        // ----------------------------------------------------------------------------------------
        // Quantification & Tximport
        // ----------------------------------------------------------------------------------------
        QUANT(reads_for_quant, use_idx)
        quant_for_multiqc = QUANT.out.quant
        quant_versions = QUANT.out.versions

        if(!params.skip_tximport){

            TXIMPORT(QUANT.out.quant.collect(), use_tx2gene)
            tximport_versions = TXIMPORT.out.versions

        } else {

            tximport_versions = Channel.empty()

        }

        // ----------------------------------------------------------------------------------------
        // MultiQC
        // ----------------------------------------------------------------------------------------

        if(!params.skip_multiqc) { MULTIQC(fastqc_for_multiqc.concat(quant_for_multiqc).collect()) }
       
    }

    // ----------------------------------------------------------------------------------------
    // Command lines and software versions
    // ----------------------------------------------------------------------------------------
    x_commands = idx_versions.concat(cat_versions, trim_versions, fastqc_versions, quant_versions, tximport_versions)
                 .map {it [1]}.flatten().collect()

    x_versions = idx_versions
                 .concat(cat_versions.first(), 
                         fastqc_versions.first(), 
                         trim_versions.first(),
                         quant_versions.first(),
                         tximport_versions)
                .map {it [0]}
                .flatten()
                .collect()

    CommandLines(x_commands, x_versions)

}

workflow { RNASEQ_PREPROCESS() }

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