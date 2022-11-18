process CommandLines {

    cpus 1

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "salmon=1.9.0" }
    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    path(cls, stageAs: "?/*")
    path(vers, stageAs: "?/*")
    
    output:
    path("commandlines.txt")
    path("versions.txt")
        
    script:
    """
    cat $cls | awk NF | grep -vE '.command.sh|> versions.txt' | sort -u > commandlines.txt
    cat $vers | awk NF | sort -u > versions.txt
    """

} 

