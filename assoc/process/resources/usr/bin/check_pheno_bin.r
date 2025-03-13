#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--listpheno"), type="character", default='init',
              help="dataset file name", metavar="character"),
  make_option(c("--listpheno_type"), type="character", default="Pheno",
              help="dataset file name", metavar="character"),
  make_option(c("--out"), type="integer", default=1,
              help="dataset file name", metavar="integer"),
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Data<-read.table(opt[['data']],header=T)
listpheno<-strsplit(opt[['listpheno']],sep=',')[[1]]
listpheno_type<-strsplit(opt[['listpheno_type']],sep=',')[[1]]
if(length(listpheno_type)==1) listpheno_type<-rep(listpheno_type,length(listpheno))

for(cmt in 1:length(listpheno)){
 if(listpheno_type[cmt]==1){
	x<-Data[,listphenop[cmt]]
	if(length(unique(x[!is.na(x)]))!=2){
            cat('pheno', listphenop[cmt], 'not binary','\n')
	    q('no',3)
	}
	minx=min(x,na.rm=T)
	maxx=max(x,na.rm=T)
	x[!is.na(x) & x==minx]<-1
	x[!is.na(x) & x==maxx]<-1
	Data[,listphenop[cmt]]<-x
 }
}
writeLines(paste(listpheno[listpheno_type==1],collapse=','), con='bin_pheno')
write.table(Data, row.names=F, col.names=T, sep='\t', quote=F, file=opt[['out']])


