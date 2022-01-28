process Quant {

    tag "$sample_id"

    label 'process_quant'

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "salmon:1.6.0" }
    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    tuple val(sample_id), path(r1), path(r2), val(libtype)
    path(idx)     
    path(tx2gene)                    

    output:
    tuple val(sample_id), path(sample_id), emit: quants
    path("$sample_id/tx2gene.txt"), emit: tx2gene
    
    script:
    def lib_type = libtype.toString().replaceAll('\\[|\\]|\'', '')
    
    // hacking to make both single and paired input work
    def reads = r2.toString() == "null" ? "-r $r1" : "-1 $r1 -2 $r2"

    """
    salmon quant --no-version-check --validateMappings \
        -i $idx -o $sample_id -l $lib_type -p $task.cpus $params.additional $reads
    
    cat $tx2gene > ${sample_id}/tx2gene.txt
    """

} 

