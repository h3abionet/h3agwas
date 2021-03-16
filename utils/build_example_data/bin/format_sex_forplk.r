#!/usr/bin/env Rscript
library("optparse")


option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-i", "--info"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-e", "--error"), type="character",
              help="dataset file name", metavar="character", default="0"),
    make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

allline<-readLines(basename(opt[['file']]))
allline<-gsub('\t$','',gsub('\t\t', '', allline,perl=T),perl=T)
DataInfo<-read.table(text=allline, sep='\t', header=T)
head(DataInfo)
Info=opt[['info']]
print(Info)
Info=strsplit(Info, split=',')[[1]]
SexHead=Info[1]
Sex1=Info[2];Sex1P=strsplit(Sex1, split=':')[[1]][1];Sex1I=strsplit(Sex1, split=':')[[1]][2]
Sex2=Info[3];Sex2P=strsplit(Sex2, split=':')[[1]][1];Sex2I=strsplit(Sex2, split=':')[[1]][2]
print(Info[4])
IID=strsplit(Info[4],split=':')[[1]][2]
errorrate<-as.numeric(opt[['error']])
nbind<-nrow(DataInfo)
nbinderror<-round(nbind*errorrate)
listind<-c()
if(errorrate>0){
print(nbinderror)
print(errorrate)
print(IID)
listind<-sample(DataInfo[,IID],nbinderror)
}


baliseerror<-DataInfo[,IID] %in% listind
DataInfo[,SexHead]<-as.character(DataInfo[,SexHead])
DataInfo$sex_plinkf<-NA
DataInfo$sex_plinkf[ !baliseerror & DataInfo[,SexHead]==Sex1I]<-Sex1P
DataInfo$sex_plinkf[!baliseerror & DataInfo[,SexHead]==Sex2I]<-Sex2P
DataInfo$sex_plinkf[baliseerror & DataInfo[,SexHead]==Sex1I]<-Sex2P
DataInfo$sex_plinkf[baliseerror & DataInfo[,SexHead]==Sex2I]<-Sex1P
write.csv(DataInfo, file=paste(opt[["out"]],'sex_change_plk.ind',sep=''))
write.table(DataInfo[,c(IID,IID, "sex_plinkf")],sep='\t', file=opt[['out']], row.names=F, col.names=F, quote=F)
