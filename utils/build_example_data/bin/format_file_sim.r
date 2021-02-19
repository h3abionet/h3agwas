#!/usr/bin/env Rscript
library("optparse")



option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Data<-read.table(opt[['file']],header=F)
for(cmt in 3:ncol(Data))Data[Data[,cmt]==-9 ,cmt]<-NA
names(Data)<-c('FID','IID', paste("pheno_",1:(ncol(Data)-2),sep=''))
write.table(Data, row.names=F, col.names=T, sep='\t', quote=F, file=opt[['out']])

