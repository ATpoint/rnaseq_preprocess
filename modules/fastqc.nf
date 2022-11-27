process FastQC {

    tag "$meta.id"

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.yml") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    tuple val(meta), path(reads, stageAs: "?/*")
            
    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script: 
    if(!meta.single_end){

        """
        fastqc -o ./ -q ${reads[0]}
        fastqc -o ./ -q ${reads[1]}

        cat .command.sh > command_lines.txt

        echo 'FastQC:' \$(fastqc --version | cut -d " " -f2) > versions.txt
        """

    } else {

        """
        fastqc -o ./ -q $reads

        cat .command.sh > command_lines.txt

        echo 'FastQC:' \$(fastqc --version | cut -d " " -f2) > versions.txt
        
        """

    }

}