#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="file contains phenotype", metavar="character"),
  make_option(c( "--ind_vcf"), type="character", default=NULL,
              help="individual list from vcf", metavar="character"),
  make_option(c( "--ind_bgen"), type="character", default=NULL,
              help="individual from bgen", metavar="character"),
  make_option(c( "--pheno"), type="character", default=NULL,
              help="phenotype separated by , ", metavar="character"),
  make_option(c( "--pheno_bin"), type="integer", default=0,
              help="phenotype is binary"),
  make_option(c("--out"), type="character",
              help="list of genes ", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

listpheno=strsplit(opt[['pheno']],split=',')[[1]]

Data<-read.table(opt[['data']], header=T)
if(any(!(listpheno %in% names(Data)))){
cat('pheno ',paste(listpheno[!(listpheno %in% names(Data))], collapse=','), 'not found in ',opt[['pheno']],'\n')
q(2,'n')
}
if(opt[['pheno_bin']]==1){
 for(pheno in listpheno){
   tb<-table(Data[,pheno])
   if(length(tb)!=2) {
    cat(tb,'\n')
    cat('pheno ',pheno, 'not binary\n')
    q(2,'n')
   }
   mintb<-min(Data[,pheno], na.rm=T)
   maxtb<-max(Data[,pheno], na.rm=T)
   if(mintb!=0 | maxtb!=1){
    newpheno<-Data[,pheno]
    newpheno[!is.na(newpheno) & newpheno==mintb]<-0
    newpheno[!is.na(newpheno) & newpheno==maxtb]<-1
    cat(mintb, maxtb,'\n')
    print(table(newpheno==mintb))
    Data[,pheno]<-newpheno
   }
 }
}
ndata<-names(Data)
Data$FIDIID<-paste(Data[,1], Data[,2],sep='_')
if(!is.null(opt[['ind_vcf']])){
 listind<-read.table(opt[['ind_vcf']])
 if(any(listind$V1 %in% listind[,1]) & any(listind$V2 %in% listind[,1])){
  write.table(Data[,ndata], file=opt[['out']], row.names=F, col.names=T, sep='\t', quote=F)
  write.table(Data[,c(1,2,1,2)], file=paste(opt[['out']],'_updateid',sep=''),  col.names=F, sep='\t', quote=F)
 }else if(any(listind$V1 %in% Data$FIDIID)){
  write.table(Data[,c('FID','IID','FIDIID','FIDIID')], file=paste(opt[['out']],'_updateid',sep=''),  col.names=F, sep='\t', quote=F, row.names=F)
  Data$FID<-Data$FIDIID
  Data$IID<-Data$FIDIID
  write.table(Data[,ndata], file=opt[['out']], row.names=F, col.names=T, sep='\t', quote=F)
}else if(any(listind$V1 %in% listind[,1])){
  write.table(Data[,c('FID','FID','FID','FID')], file=paste(opt[['out']],'_updateid',sep=''),  col.names=F, sep='\t', quote=F, row.names=F)
  Data$IID<-Data$FID
  write.table(Data[,ndata], file=opt[['out']], row.names=F, col.names=T, sep='\t', quote=F)
 }else{
  cat(Data$FID)
  cat('\n')
  cat(Data$IID)
  cat('no same individual between vcf and pheno file\nexit\n')
  q('n', 2)
 }
}else if(!is.null(opt[['ind_bgen']])){
  DataSample<-read.table(opt[['ind_bgen']], header=T)
  DataSample<-Data[-1,]
  DataTmp<-merge(Data, DataSample)
  if(nrow(DataTmp)>1){
     write.table(Data[,ndata], file=opt[['out']], row.names=F, col.names=T, sep='\t', quote=F)
     write.table(Data[,c(1,2,1,2)], file=paste(opt[['out']],'_updateid',sep=''),  col.names=F, sep='\t', quote=F)
  }else if(any(DataTmp[,1] %in%  Data$FIDIID)){
   write.table(Data[,c('FID','IID','FIDIID','FIDIID')], file=paste(opt[['out']],'_updateid',sep=''),  col.names=F, sep='\t', quote=F, row.names=F)
   Data$FID<-Data$FIDIID
   Data$IID<-Data$FIDIID
   write.table(Data[,ndata], file=opt[['out']], row.names=F, col.names=T, sep='\t', quote=F)
 } else{
  cat(Data$FID)
  cat('\n')
  cat(Data$IID)
  cat('no same individual between vcf and pheno file\nexit\n')
  q('n', 2)
 }

}else{
     write.table(Data[,ndata], file=opt[['out']], row.names=F, col.names=T, sep='\t', quote=F) 
     write.table(Data[,c(1,2,1,2)], file=paste(opt[['out']],'_updateid',sep=''),  col.names=F, sep='\t', quote=F)
}

