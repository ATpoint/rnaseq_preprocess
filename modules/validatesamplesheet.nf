process ValidateSamplesheet {

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    path(samplesheet)
            
    output:
    path("samplesheet_validated.csv"), emit: samplesheet
    
    script: 
    """
    Rscript --vanilla $baseDir/bin/validate_samplesheet.R $samplesheet
    """     

}