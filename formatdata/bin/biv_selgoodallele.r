#!/usr/bin/env Rscript
library(data.table)
args = commandArgs(trailingOnly=TRUE)
databim<-fread(args[1])
databim$chrpos<-paste(databim$V1, databim$V4)
tmp=table(databim$chrpos)
tmp2<-table(as.character(databim$V2))
balise<-(as.character(databim$V5) %in% c("A", "T", "C", "G")) & (as.character(databim$V6) %in% c("A", "T", "C", "G"))  & (databim$chrpos %in% names(tmp[tmp==1])) & (databim$V2 %in% tmp2[tmp2==1])
writeLines(as.character(databim$V2[!balise]), con=args[2])

