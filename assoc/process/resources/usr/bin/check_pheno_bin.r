#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--listpheno"), type="character", default='init',
              help="dataset file name", metavar="character"),
  make_option(c("--listpheno_bin"), type="character", default="Pheno",
              help="dataset file name", metavar="character"),
  make_option(c("--software"), type="character", default="plink",
              help="dataset file name", metavar="character"),
  make_option(c("--out"), type="character", default="out",
              help="dataset file name", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Data<-read.table(opt[['data']],header=T)
listpheno<-strsplit(opt[['listpheno']],split=',')[[1]]
listpheno_type<-strsplit(opt[['listpheno_bin']],split=',')[[1]]
if(length(listpheno_type)==1) listpheno_type<-rep(listpheno_type,length(listpheno))
if(opt[['software']]=='plink'){
minsoftware=1
maxsoftware=2
}else {
minsoftware=0
maxsoftware=1
}

for(cmt in 1:length(listpheno)){
 if(listpheno_type[cmt]==1){
	x<-Data[,listpheno[cmt]]
	if(length(unique(x[!is.na(x)]))!=2){
            cat('pheno', listpheno[cmt], 'not binary','\n')
	    q('no',3)
	}
	minx=min(x,na.rm=T)
	maxx=max(x,na.rm=T)
	x[!is.na(x) & x==minx]<-minsoftware
	x[!is.na(x) & x==maxx]<-maxsoftware
	Data[,listpheno[cmt]]<-x
 }
}
writeLines(paste(listpheno[listpheno_type==1],collapse=','), con='bin_pheno')
write.table(Data, row.names=F, col.names=T, sep='\t', quote=F, file=opt[['out']])


