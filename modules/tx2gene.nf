process Tx2Gene {

    label 'process_tx2gene'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    path(gtf)
            
    output:
    path("tx2gene.txt"), emit: tx2gene
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
        
    script:    
    """
    Rscript --vanilla ${baseDir}/bin/tx2gene.R $gtf tx2gene.txt $params.transcript_id $params.transcript_name $params.gene_id $params.gene_name $params.gene_type

    cat .command.sh > command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    echo 'rtracklayer:' \$(Rscript -e "cat(as.character(packageVersion('rtracklayer')))") >> versions.txt
    """                

}