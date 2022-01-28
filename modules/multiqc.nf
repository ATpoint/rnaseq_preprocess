process MultiQC {

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda params.environment }
    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

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