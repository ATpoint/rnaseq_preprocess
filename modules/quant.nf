// Quantify fastq files with salmon

process SalmonQuant {

    tag "$sample_id"

    cpus   params.threads
    memory params.mem

    publishDir params.outdir, mode: params.pubmode

    input:
    tuple val(sample_id), path(reads)
    path(idx)                         

    output:
    tuple val(sample_id), path(sample_id), emit: quants
    
    script:

    def readfiles = (params.mode=='single') ? "-r $reads" : "-1 ${reads[0]} -2 ${reads[1]}"

    """
    salmon quant --no-version-check --validateMappings \
        -i $idx -o $sample_id -l $params.libtype -p $task.cpus $params.additional $readfiles
    """

} 

