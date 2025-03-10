#!/usr/bin/env Rscript
library("optparse")
#     addpheno_pcs.r --data $phenofile --pcs $head".eigenvec" --out newfilepheno
option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--pcs"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--covar"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Data<-read.table(opt[['data']] ,header=T)
Pcs<-read.table(opt[['pcs']])
names(Pcs)<-c('FID', 'IID', paste('Pcs_',1:(ncol(Pcs)-2),sep=''))
DataPcs<-merge(Data, Pcs, by=c(1,2))
write.table(DataPcs, file=opt[['out']], sep='\t', col.names=T, row.names=F, quote=F)

listcovar=""
if(!is.null(opt[['covar']]))listcovar=paste(opt[['covar']],",",sep='')
listpcs<-paste(paste('Pcs_',1:(ncol(Pcs)-2),sep=''), collapse=',')

writeLines(paste(listcovar,listpcs,sep=''), con=paste(opt[['out']], '.covar', sep=''))



