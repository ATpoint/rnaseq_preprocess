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
    tuple val(sampleid), path(r1, stageAs: "?/*"), path(r2, stageAs: "?/*")
            
    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    
    script: 
    // hacking to make both single and paired input work
    def reads = r2.baseName.toString() == "null" ? r1 : "$r1 $r2"
    """
    fastqc --threads 1 -o ./ $reads
    """     

}