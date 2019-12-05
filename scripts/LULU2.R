setwd(".")

require(lulu)
require(methods)
require(dplyr)

args <- commandArgs()
dist <- args[6]
tsv <- args[7]
thr <- args[8]

curtsv <- gsub(".tsv", "_curated.shared", tsv)

matchlist = read.delim(dist, header=FALSE,as.is=TRUE, stringsAsFactors=FALSE)
otus = read.delim(tsv,row.names=1,header=T,sep="\t")
curated_result <- lulu(as.data.frame(t(otus)), matchlist, minimum_match = thr)
lulus = t(curated_result$curated_table)
Group <- rownames(lulus)
numOtus <- ncol(lulus)
rownames(lulus) <- NULL
data <- cbind(label="swarm_1", Group,numOtus,lulus)

write.table(data,curtsv, row.names=FALSE, quote=F, sep="\t")
