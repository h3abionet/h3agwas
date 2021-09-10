#!/usr/bin/env Rscript
library("optparse")



option_list = list(
  make_option(c("--metal"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--resume"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--out"), type="character", default=NULL,
              help="dataset file name", metavar="character")

)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


DataMetal<-read.table(opt[['metal']],header=T)
DataResume<-read.csv(opt[['resume']])
DataResume$MarkerName<-paste(DataResume$CHR,DataResume$BP,sep='_')
DataM<-merge(DataResume,DataMetal, by='MarkerName',all.x=T)
write.csv(DataM, row.names=F, file=opt[['out']])

