process Quant {

    tag "$sample_id"

    label 'process_quant'

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") ? null : filename } 
    ]

    if(workflow.profile.contains('conda'))  { conda "salmon=1.9.0" }
    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    tuple val(sample_id), path(r1, stageAs: "?/*"), path(r2, stageAs: "?/*"), val(libtype)
    path(idx)     
    path(tx2gene)                    

    output:
    path(sample_id), emit: quants
    path("$sample_id/tx2gene.txt.gz"), emit: tx2gene
    path("commandlines.txt"), emit: commandlines
    path("versions.txt"), emit: versions
    
    script:

    def lib_type = libtype[0]
    
    // hacking to make both single and paired input work
    def reads = r2[0].baseName.toString().matches("null") ? "-r $r1" : "-1 $r1 -2 $r2"

    // Make sure gcBias not used in single-end data
    def additional = r2[0].baseName.toString().matches("null") ? params.quant_additional.replaceAll('--gcBias', '') : params.quant_additional

    """
    salmon quant --no-version-check --validateMappings -i $idx -o $sample_id -l $lib_type -p $task.cpus $additional $reads
    cat $tx2gene | gzip > ${sample_id}/tx2gene.txt.gz
    cat .command.sh | awk NF | grep -v '^#!' > commandlines.txt
    echo 'salmon:' \$(salmon --version | cut -d " " -f2) > versions.txt
    """

}