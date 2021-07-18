#/ Prepare a tx2gene map for salmon/tximport:

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 7) { 
  stop("Usage: tx2gene.R <gtf> <outname> <tx_id> <tx_name> <gene_id> <gene_name> <gene_type>",
       call.=FALSE)
}
gtf <- args[1]; outname <- args[2]; transcript_id <- args[3]; transcript_name <- args[4]
gene_id <- args[5]; gene_name <- args[6]; gene_type <- args[7]

tx2gene <- data.frame(rtracklayer::import(gtf))
tx2gene <- unique(tx2gene[tx2gene$type=="transcript",])
tx2gene <- tx2gene[,c(transcript_id, gene_id, gene_name, gene_type, transcript_name)]

tx2gene <- data.frame(transcript_id=tx2gene[,transcript_id],
                      name=paste(tx2gene[,gene_id], tx2gene[,gene_name], sep="_"),
                      tx2gene[,2:ncol(tx2gene)])

write.table(x=tx2gene, file=outname, sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE)
