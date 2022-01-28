process MultiQC {

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "multiqc=1.11"}
    if(workflow.profile.contains('docker')) { container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0" }
    if(workflow.profile.contains('singularity')) { container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0" }

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