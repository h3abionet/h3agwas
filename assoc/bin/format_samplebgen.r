#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--sample"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c( "--out_sample"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c( "--out_samplesaige"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c( "--data"), type="character", default=NULL,
              help="dataset file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

DataPheno<-read.table(opt[['data']], header=T)
DataSample<-read.table(opt[['sample']], header=T)
DataSample[,1]<-as.character(DataSample[,1]);DataSample[,2]<-as.character(DataSample[,2])
DataPheno[,1]<-as.character(DataPheno[,1]);DataPheno[,2]<-as.character(DataPheno[,2])
balise<-which((DataSample[,1] %in% DataPheno[,1]) & (DataSample[,2] %in% DataPheno[,2]))

if(length(balise)<5){
DataSample[,1]<-gsub("^[A-Za-z0-9]+_", "", DataSample[,1])
DataSample[,2]<-gsub("_[A-Za-z0-9]+$", "", DataSample[,2])
}
balise<-which((DataSample[,1] %in% DataPheno[,1]) & (DataSample[,2] %in% DataPheno[,2]))
if(length(balise)<5){
cat("no found common individual between data and bgen")
q(2)
}

DataSampleSaige<-DataSample
DataSampleSaige$sex<-c("D", rep(0,nrow(DataSampleSaige)-1))

write.table(DataSample, file=opt[['out_sample']], row.names=F, col.names=T, quote=F)
write.table(DataSampleSaige, file=opt[['out_samplesaige']], row.names=F, col.names=T, quote=F)


