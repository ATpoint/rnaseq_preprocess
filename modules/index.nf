process Idx {

    label 'process_idx'

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "salmon:1.6.0"}
    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

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
    gzip -cd $genome | grep '^>' | cut -d " " -f 1 | awk '{gsub(\">\",\"\");print}' > $decoynames

    cat $txtome $genome > $gentrome

    salmon index --no-version-check \
        -t $gentrome \
        -d $decoynames \
        -i $idxname \
        -p $task.cpus \
        $params.additional
    """                

}