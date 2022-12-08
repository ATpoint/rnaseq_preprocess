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

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'MultiQC:' \$(multiqc --version 2>&1 | cut -d " " -f3)  > versions.txt
    """     

}