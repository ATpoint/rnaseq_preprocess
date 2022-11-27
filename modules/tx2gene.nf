process Tx2Gene {

    label 'process_tx2gene'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    path(gtf)
            
    output:
    path("tx2gene.txt"), emit: tx2gene
        
    script:    
    """
    Rscript --vanilla ${baseDir}/bin/tx2gene.R \
        $gtf tx2gene.txt $params.transcript_id $params.transcript_name $params.gene_id $params.gene_name $params.gene_type
    """                

}