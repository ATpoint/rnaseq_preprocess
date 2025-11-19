process Tximport {

    label 'process_tximport'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container
    
    input:
    path(quants) 
    path(tx2gene)
        
    output:
    path("counts.txt.gz")
    path("lengths.txt.gz")
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
                
    script: 
    def q = quants.join(',').toString()
    """
    Rscript --vanilla $baseDir/bin/tximport.R $q $tx2gene

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    echo 'tximport:' \$(Rscript -e "cat(as.character(packageVersion('tximport')))") >> versions.txt
    """      

}
