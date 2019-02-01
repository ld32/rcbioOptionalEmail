#!/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2) {
  stop("Usage: change.tmp.table.columns.to.match.rsem.result.r  [tpm.file] [rsem.result]")
}


if(!file.exists(args[1])){
   stop(paste("File not exit: ", args[1]))
}

if(!file.exists(args[2])){
   stop(paste("File not exit: ", args[2])) 
}

tmp=read.table(args[1] ,check.name=FALSE)

rsem=read.table(args[2],check.name=FALSE)

new.tmp=tmp[,colnames(rsem)]

write.table(new.tmp, file=paste("correct.header.", args[1], sep=""))

