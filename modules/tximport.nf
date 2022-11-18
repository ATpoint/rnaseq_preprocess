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
    path("versions.txt"), emit: versions
                
    script: 
    def q = quants.join(',').toString()
    """
    Rscript --vanilla $baseDir/bin/tximport.R $q $tx2gene
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    echo 'tximport:' \$(Rscript -e "cat(as.character(packageVersion('tximport')))") >> versions.txt
    """      

}