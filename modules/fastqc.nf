process FastQC {

    tag "$sampleid"

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda params.environment }
    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    tuple val(sampleid), path(r1, stageAs: "?/*"), path(r2, stageAs: "?/*"), val(libtype)
            
    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    
    script: 
    r1_use = r1
    r2_use = r2.baseName.toString().contains("null") ? 'false' : 'true'
    
    """
    fastqc --threads 1 -o ./ $r1
    if [[ $r2_use == 'true' ]]; then fastqc --threads 1 -o ./ $r2; fi
    """     

}