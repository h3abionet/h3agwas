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
  make_option(c("--covar_type"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--npcs"), type="integer", default=NULL,
              help="dataset file name"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Data<-read.table(opt[['data']] ,header=T)
Pcs<-read.table(opt[['pcs']]);Pcs<-Pcs[,1:(opt[['npcs']]+2)]
names(Pcs)<-c('FID', 'IID', paste('Pcs_',1:(ncol(Pcs)-2),sep=''))
DataPcs<-merge(Data, Pcs, by=c(1,2))
write.table(DataPcs, file=opt[['out']], sep='\t', col.names=T, row.names=F, quote=F)

listcovar=""
if(!is.null(opt[['covar']]))listcovar=paste(opt[['covar']],",",sep='')
listpcs<-paste(paste('Pcs_',1:(ncol(Pcs)-2),sep=''), collapse=',')

writeLines(paste(listcovar,listpcs,sep=''), con=paste(opt[['out']], '.covar', sep=''))

covar_n<-strsplit(opt[['covar']],split=',')[[1]]
covar_n<-length(covar_n[nchar(covar_n)>0])

listcovar_pcs<-rep(0,opt[['npcs']])
listcovar_type<-c()
if(!is.null(opt[['covar']]))listcovar_type<-c('0')
if(!is.null(opt[['covar_type']]))listcovar_type<-strsplit(opt[['covar_type']],split=',')[[1]]
if(!is.null(opt[['covar']]) & length(listcovar_type)==1){
	listcovar_type<-rep(listcovar_type, covar_n)
}
if(length(listcovar_type)!=covar_n){
cat('length(listcovar_type)!=covar_n') 
q(5)
}
listcovar_type<-c(listcovar_type,listcovar_pcs)
writeLines(paste(listcovar_type,collapse=','), con=paste(opt[['out']], '.covartype', sep=''))







