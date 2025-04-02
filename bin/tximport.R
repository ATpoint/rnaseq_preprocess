library(tximport)

args = commandArgs(trailingOnly=TRUE)
if(length(args) != 2) { 
  stop("Usage: tximport.R <salmon_dir> <tx2gene>",call.=FALSE)
}

quants   <- args[1]
quants   <- strsplit(quants, split=",")[[1]]
bname    <- quants
tx2gene  <- read.delim(args[2], header=TRUE)
q        <- paste0(quants, "/quant.sf.gz")
names(q) <- bname
q        <- q[order(names(q))]

# Add a simple check if there is a mismatch between quant.sf and tx2gene that can be solved
# be removing txversion or bar
q1 <- read.delim(q[1], header=TRUE)
v <- q1[,1,drop=TRUE]
as.is <- v %in% tx2gene[,1,drop=TRUE] 
no.dot <- gsub("\\..*", "", v) %in% tx2gene[,1,drop=TRUE] 
no.bar <- gsub("\\|.*", "", v) %in% tx2gene[,1,drop=TRUE] 
l <- c(sum(as.is), sum(no.dot), sum(no.bar))
w <- which(l==max(l))

if(1 %in% w){
  ignoreTxVersion <- FALSE
  ignoreAfterBar  <- FALSE
}

if(2 %in% w & length(w)==1){
  ignoreTxVersion <- TRUE
  ignoreAfterBar  <- FALSE
}

if(3 %in% w & length(w)==1){
  ignoreTxVersion <- FALSE
  ignoreAfterBar  <- TRUE
}

txi <- tximport::tximport(files=q, tx2gene=tx2gene, 
                          type="salmon", importer=read.delim,
                          countsFromAbundance="lengthScaledTPM",
                          dropInfReps=TRUE, 
                          ignoreTxVersion=ignoreTxVersion, 
                          ignoreAfterBar=ignoreAfterBar)

#/ save counts with length being the median of average tx length:
tgz <- gzfile("counts.txt.gz", "w")
write.table(x=data.frame(gene=rownames(txi$counts), txi$counts), 
            file=tgz, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
close(tgz)

tgz <- gzfile("lengths.txt.gz", "w")
write.table(x=data.frame(gene=rownames(txi$length), txi$counts), 
            file=tgz, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
close(tgz)

tgz <- gzfile("tx2gene.txt.gz", "w")
write.table(x=tx2gene, 
            file=tgz, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
close(tgz)


