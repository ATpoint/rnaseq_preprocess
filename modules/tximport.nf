process Tximport {

    label 'process_tximport'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda params.environment }
    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }
    
    input:
    path(quants) 
    path(tx2gene)
        
    output:
    path("counts.txt.gz")
    path("lengths.txt.gz")
                
    script: 
    def q = quants.join(',').toString()
    """
    Rscript --vanilla $baseDir/bin/tximport.R $q $tx2gene
    """      

}