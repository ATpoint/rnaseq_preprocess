process ValidateSamplesheet {

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    container params.container

    input:
    path(samplesheet)
            
    output:
    path("samplesheet_validated.csv"), emit: samplesheet
    
    script: 
    """
    Rscript --vanilla $baseDir/bin/validate_samplesheet.R $samplesheet

    echo ${task.process}: > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    """     

}