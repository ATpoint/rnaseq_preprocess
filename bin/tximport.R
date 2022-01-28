library(tximport)

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 3) { 
  stop("Usage: tximport.R <salmon_dir> <outname> <tx2gene>",call.=FALSE)
}

quants   <- args[1]
quants   <- strsplit(quants, split=",")[[1]]
bname    <- quants
outname  <- args[2]
tx2gene  <- read.delim(args[3], header=TRUE)
q        <- paste0(quants, "/quant.sf")
names(q) <- bname

txi <- tximport::tximport(files=q, tx2gene=tx2gene, 
                          type="salmon", importer=read.delim,
                          countsFromAbundance="lengthScaledTPM")

q#/ save counts with length being the median of average tx length:
tgz <- gzfile(paste0(outname, "counts_genelevel.txt.gz"), "w")

rl <- apply(txi$length, 1, median)

write.table(x=data.frame(gene=rownames(txi$counts), 
                         length=ceiling(rl),
                         txi$counts), 
            file=tgz, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

close(tgz)
                     
                     