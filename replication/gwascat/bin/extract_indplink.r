#!/usr/bin/env Rscript
library("optparse")
#     addpheno_pcs.r --data $phenofile --pcs $head".eigenvec" --out newfilepheno
option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--pheno"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--cov"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Data<-read.table(opt[['data']] ,header=T)

listdata=names(Data)[1:2]
if(!is.null(opt[['pheno']])){
listdata<-c(listdata, opt[['pheno']])
}

if(!is.null(opt[['cov']])){
listdata<-c(listdata, strsplit(opt[['cov']], split=',')[[1]])
}
Data2<-na.omit(Data[,listdata])
write.table(Data2[,1:2], file=opt[['out']], sep='\t', row.names=F, col.names=F, quote=F)
