library(tximport)

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 2) { 
  stop("Usage: tximport.R <salmon_dir> <tx2gene>",call.=FALSE)
}

quants   <- args[1]
quants   <- strsplit(quants, split=",")[[1]]
bname    <- quants
tx2gene  <- read.delim(args[2], header=TRUE)
q        <- paste0(quants, "/quant.sf")
names(q) <- bname
q        <- q[order(names(q))]

txi <- tximport::tximport(files=q, tx2gene=tx2gene, 
                          type="salmon", importer=read.delim,
                          countsFromAbundance="lengthScaledTPM",
                          dropInfReps=TRUE)

#/ save counts with length being the median of average tx length:
tgz <- gzfile("counts.txt.gz", "w")
write.table(x=data.frame(gene=rownames(txi$counts), txi$counts), 
            file=tgz, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
close(tgz)

tgz <- gzfile("lengths.txt.gz", "w")
write.table(x=data.frame(gene=rownames(txi$length), txi$counts), 
            file=tgz, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
close(tgz)
                     
                     