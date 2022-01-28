process Tximport {

    label 'process_tximport'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda params.environment }
    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }
    
    input:
    path(quants) 
    val(outname)  // collected paths to all salmon dirs
    path(tx2gene)
        
    output:
    path("*counts_genelevel.txt.gz")
    //path("*_infreps.txt.gz"), optional: true
            
    script: 
    def q = quants.join(',').toString()
    """
    Rscript --vanilla $baseDir/bin/tximport.R $q $outname $tx2gene
    """      

}