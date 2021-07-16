// Create genome-decoyed index with salmon and a tx2gene mapping table

process SalmonIndex {

    cpus   params.threads
    memory params.mem

    publishDir params.outdir, mode: params.pubmode

    input:
    path(txtome)
    path(genome) 
    val(idxname)
    path(gtf)
        
    output:
    path("decoynames.txt")
    path("gentrome.fa.gz")
    path(idxname), emit: idx
    path("tx2gene.txt"), emit: tx2gene
    
    script: 

    def decoynames  = "decoynames.txt"
    def gentrome    = "gentrome.fa.gz"

    """
    zgrep '^>' $genome | cut -d " " -f 1 | awk '{gsub(\">\",\"\");print}' > $decoynames

    cat $txtome $genome > $gentrome

    salmon index --no-version-check \
        -t $gentrome \
        -d $decoynames \
        -i $idxname \
        -p $task.cpus \
        $params.additional

    Rscript --vanilla ${baseDir}/src/tx2gene.R \
        annot.gtf.gz tx2gene.txt $params.transcript_id $params.transcript_name $params.gene_id $params.gene_name $params.gene_type

    """                

}