process Quant {

    tag "$meta.id"

    label 'process_quant'

    errorStrategy 'finish'
    
    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    container params.container

    input:
    tuple val(meta), path(reads)
    path(idx)     

    output:
    path("$meta.id"), emit: quant
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script:
    def lib_type = meta.lib_type
    
    // Make sure --gcBias not used in single-end data
    def additional = meta.single_end ? params.quant_additional.replaceAll('--gcBias', '') : params.quant_additional

    use_reads = meta.single_end ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"

    """
    salmon quant --no-version-check --validateMappings -i $idx -o ${meta.id} -l $lib_type -p $task.cpus $additional $use_reads
    gzip --best ${meta.id}/quant.sf

    echo ${task.process}:${meta.id} > command_lines.txt
    cat .command.sh | grep -vE '^#!/bin|versions.txt\$|command_lines.txt\$|cat \\.command.sh' | sed 's/  */ /g' | awk NF >> command_lines.txt

    echo 'salmon:' \$(salmon --version | cut -d " " -f2) > versions.txt

    """

}