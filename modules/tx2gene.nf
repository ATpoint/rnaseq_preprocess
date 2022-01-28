process Tx2Gene {

    label 'process_tx2gene'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "bioconductor-rtracklayer=1.54.0"}
    if(workflow.profile.contains('docker')) { container "quay.io/biocontainers/bioconductor-rtracklayer:1.54.0--r41ha2fdcc6_1" }
    if(workflow.profile.contains('singularity')) { container "quay.io/biocontainers/bioconductor-rtracklayer:1.54.0--r41ha2fdcc6_1" }

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