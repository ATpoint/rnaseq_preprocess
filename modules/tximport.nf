process Tximport {

    label 'process_tximport'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }
    
    input:
    path(quants) 
    path(tx2gene)
        
    output:
    path("counts.txt.gz")
    path("lengths.txt.gz")
    path("tx2gene.txt.gz")
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
                
    script: 
    def q = quants.join(',').toString()
    """
    Rscript --vanilla $baseDir/bin/tximport.R $q $tx2gene
    gzip -c $tx2gene > tx2gene.txt.gz

    cat .command.sh > command_lines.txt
    
    echo 'R:' \$(R --version | head -n1 | cut -d " " -f3) > versions.txt
    echo 'tximport:' \$(Rscript -e "cat(as.character(packageVersion('tximport')))") >> versions.txt
    """      

}