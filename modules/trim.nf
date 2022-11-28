process Trim {

    tag "$meta.id"

    cpus   1
    memory 1.GB

    errorStrategy 'finish'

    publishDir = [
        path: params.outdir,
        mode: params.publishmode,
        saveAs: { filename -> !(filename.endsWith("fq.gz") & params.keep) ||
                               filename.equals("versions.txt") || filename.equals("command_lines.txt") ? null : filename } 
    ]

    if(workflow.profile.contains('docker')) { container params.container }
    if(workflow.profile.contains('singularity')) { container params.container }

    input:
    tuple val(meta), path(reads, stageAs: "?/*")
            
    output:
    tuple val(meta), path("*.fq.gz"), emit: fastq_tuple
    tuple path("versions.txt"), path("command_lines.txt"), emit: versions
    
    script: 
    def compress_level = params.keep ? '-6' : '--fast'
    if(!meta.single_end){

        """
        seqtk trimfq -L $params.trim_length ${reads[0]} | gzip $compress_level > ${meta.id}_R1_trimmed.fq.gz
        seqtk trimfq -L $params.trim_length ${reads[1]} | gzip $compress_level > ${meta.id}_R2_trimmed.fq.gz

        cat .command.sh > command_lines.txt

        echo 'seqtk:' \$(seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d " " -f2) > versions.txt
        """

    } else {

        """
        seqtk trimfq -L $params.trim_length ${reads} | gzip $compress_level > ${meta.id}_R1_trimmed.fq.gz

        cat .command.sh > command_lines.txt

        echo 'seqtk:' \$(seqtk 2>&1 | head -n 3 | tail -n 1 | cut -d " " -f2) > versions.txt
        
        """

    }

}