process Tximport {

    label 'process_tximport'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "bioconductor-tximport=1.22.0"}
    if(workflow.profile.contains('docker')) { container "quay.io/biocontainers/bioconductor-tximport:1.22.0--r41hdfd78af_0" }
    if(workflow.profile.contains('singularity')) { container "quay.io/biocontainers/bioconductor-tximport:1.22.0--r41hdfd78af_0" }
    
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