process CommandLines {

    cpus 1

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    container params.container

    input:
    path(commands, stageAs: "?/*")
    path(versions, stageAs: "?/*")
    
    output:
    path("command_lines.txt")
    path("software_versions.txt")
        
    script:
    """
    cat $commands | awk NF | grep -vE '.command.sh|> versions.txt|^#!' | sort --ignore-case -u > command_lines.txt
    cat $versions | awk NF | sort --ignore-case -u > software_versions.txt
    """

} 

