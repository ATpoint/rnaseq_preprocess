process MultiQC {

    cpus   1
    memory 1.GB

    label 'process_multiqc'

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    container params.container

    input:
    path(everything)
            
    output:
    path "*multiqc_report.html"
    path "*_data" 
    
    script: 
    """
    multiqc .
    """     

}