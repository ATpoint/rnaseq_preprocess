#/ aggregate for the transcript-level abundance estimates to the gene level

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) { 
  stop("Usage: tximport.R <salmon_dir> <outname> <tx2gene>",call.=FALSE)
}

quants <- args[1]
outname <- args[2]
tx2gene <- args[3]

tx2gene <- read.delim(tx2gene, header=TRUE)

q <- paste0(quants, "/quant.sf")
names(q) <- outname

txi <- tximport::tximport(files=q, tx2gene=tx2gene, type="salmon",
                          importer=read.delim)

write.table(data.frame(Gene=rownames(txi$counts), txi$counts),
            paste0(outname, "_counts.txt"), sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE)

write.table(data.frame(Gene=rownames(txi$length), txi$length),
            paste0(outname, "_lengths.txt"), sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE)

if("infReps" %in% names(txi)){
  txi$infReps <- txi$infReps[[1]]
  fgz <- gzfile(paste0(outname, "_infreps.txt.gz"), "w")
  write.table(data.frame(Gene=rownames(txi$infReps), txi$infReps),
              file=fgz, sep="\t",
              col.names=TRUE, row.names=FALSE, quote=FALSE)
  close(fgz)
}
                     
                     