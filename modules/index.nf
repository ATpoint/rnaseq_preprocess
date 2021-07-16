// Create genome-decoyed index with salmon 

process SalmonIndex {

    cpus   params.threads
    memory params.mem

    publishDir params.outdir, mode: params.pubmode

    input:
    path(txtome)
    path(genome) 
    val(idxname)
        
    output:
    path("decoynames.txt")
    path("gentrome.fa.gz")
    path(idxname), emit: idx
    
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
    """                

}