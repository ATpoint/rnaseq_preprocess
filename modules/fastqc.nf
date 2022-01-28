process FastQC {

    tag "$sampleid"

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir params.outdir, mode: params.publishmode

    if(workflow.profile.contains('conda'))  { conda "fastqc=0.11.9"}
    if(workflow.profile.contains('docker')) { container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1" }
    if(workflow.profile.contains('singularity')) { container "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1" }


    input:
    tuple val(sampleid), path(r1), path(r2)
            
    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    
    script: 
    // hacking to make both single and paired input work
    def reads = r2.toString() == "null" ? r1 : "$r1 $r2"
    """
    fastqc --threads 1 $reads
    """     

}