#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--bfile"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c( "--chro_cond"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--pos_cond"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--data"), type="character", default="out.txt",
              help="out prefix", metavar="character"),
  make_option(c("--out"), type="character",
              help="list of genes ", metavar="character")
);

extract_pos<-function(bfile,listrs, tmpf='tmp'){

writeLines(listrs, con=paste(tmpf,'.rs', sep=''))
syst<-paste('plink -bfile ', bfile, ' --extract ', tmpf, '.rs --recode -out ', tmpf, sep='')
system(syst)
Data<-read.table(paste(tmpf, '.ped',sep=''),sep=' ')
DataMap<-read.table(paste(tmpf, '.map',sep=''))
Cmt2<-1
Data2<-Data[,c(1,2)]
listrsname<-c()
for(Cmt in seq(7, ncol(Data),2)){
rsname<-as.character(DataMap$V2[Cmt2])
lref<-c(Data[,Cmt],Data[,Cmt+1]);lref<-unique(lref[lref!=-9]);ref<-lref[1];alt<-lref[2]
Data2[,rsname]<-NA
Data2[Data[,Cmt]==Data[,Cmt+1] & Data[,Cmt]==ref,rsname]<-0
Data2[Data[,Cmt]==Data[,Cmt+1] & Data[,Cmt]==alt,rsname]<-1
Data2[Data[,Cmt]!=Data[,Cmt+1],rsname]<-0.5
Cmt2<-Cmt2+1
}
return(Data2)
}
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tmpf='tmp'
bfile=opt[['bfile']]
bimfile=paste(bfile,'.bim',sep='')
listpos=as.integer(strsplit(opt[['pos_cond']],split=',')[[1]])
DataBim<-read.table(bimfile)
listrs=as.character(DataBim$V2[DataBim$V4 %in% listpos])

filepheno<-opt[['data']]
DataI<-extract_pos(bfile,listrs)
DataPheno<-read.table(filepheno, header=T)
DataPheno2<-merge(DataPheno, DataI, by=c(1,2),al.x=T)
write.table(DataPheno2, file=opt[['out']],sep='\t', row.names=F, col.names=T, quote=F)
