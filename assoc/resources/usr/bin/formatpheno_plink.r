#!/usr/bin/env Rscript
library(ggplot2)
library("optparse")

option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--balise"), type="character", default='init',
              help="dataset file name", metavar="character"),
  make_option(c("--pheno"), type="character", default="Pheno",
              help="dataset file name", metavar="character"),
  make_option(c("--covar"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--binary"), type="integer", default=1,
              help="dataset file name", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


Data<-read.table(opt[['data']], header=T)
headbalise<-opt[['balise']]
pheno<-opt[['pheno']]
covar<-opt[['covar']]
if(!is.null(covar)){
covar<-strsplit(covar, split=',')[[1]]
}
if(headbalise=='init') DataTmp<-na.omit(Data[,c(names(Data)[1:2], pheno,covar)]) else{
 Data<-na.omit(Data[,c(names(Data)[1:2], pheno,covar, headbalise)])
 DataTmp<-Data[Data[,headbalise]==T,]
}
pheno2<-pheno
if(opt[["binary"]]==1){
unix<-length(as.integer(unique(DataTmp[,'Pheno'])))
print(length(unix))
if(unix!=2){
print('more than 2 pheno ')
cat(unix)
q(save = "no", status = 2)
}
uni<-unique(DataTmp[,'Pheno'])
DataTmp$Pheno2[DataTmp$Pheno==min(uni)]<-1
DataTmp$Pheno2[DataTmp$Pheno==max(uni)]<-2
pheno2="Pheno2"
}
DataTmp<-DataTmp[,c('FID', 'IID',pheno2,covar)];names(DataTmp)<-c('FID', 'IID',pheno,covar)
write.table(DataTmp[,c('FID', 'IID',pheno,covar)], file=opt[['out']],sep='\t',quote=F, row.names=F, col.names=T)
