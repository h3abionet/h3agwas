#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--ld"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--list_file_assoc"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--file_I"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--file_merge"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--pos_ref"), type="integer", default=NULL,
              help="dataset file name", metavar="character"),
    make_option(c("--out"), type="character", default="out.svg",
              help="output file name [default= %default]", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


fileld=opt[['ld']]
posref=opt[['pos_ref']]

DataLD<-read.table(fileld, header=T)
DataLD2<-DataLD[,c('CHR_B','BP_B','SNP_B','CHR_A','BP_A','SNP_A','R2')];names(DataLD2)<-names(DataLD)
DataLDAll<-rbind(DataLD,DataLD2)

DataLDAll<-DataLDAll[DataLDAll$BP_B==posref,c('CHR_A','BP_A','SNP_A','R2')]


listfile=strsplit(opt[['list_file_assoc']], split=',')[[1]]
Cmt<-1
for(File in listfile){
headcond<-gsub('_',':',gsub('.assoc.txt','',basename(File)))
Data<-read.table(File , header=T)
Data$rscond<-headcond
if(Cmt==1)DataEachRs<-Data
else DataEachRs<-rbind(DataEachRs,Data)
Cmt<-Cmt+1
}

DataI<-read.table(opt[['file_I']] ,header=T);DataI$rscond<-'Init'
DataMerge<-read.table(opt[['file_merge']], header=T);DataMerge$rscond<-'Merge'

DataEachRs<-rbind(DataI,DataEachRs,DataMerge)

DataEachRs_posref<-DataEachRs[DataEachRs$ps==posref,]
DataEachRs_posref$Num<-1:nrow(DataEachRs_posref)

print(DataEachRs_posref$rscond)
print(DataLDAll$SNP_A)
DataEachRs_posref<-merge(DataEachRs_posref ,DataLDAll, by.x='rscond', by.y='SNP_A', all.x=T)
DataEachRs_posref<-DataEachRs_posref[order(DataEachRs_posref$Num),]
DataEachRs_posref$R2[1]<-1
pdf(paste(opt[['out']],'.pdf',sep=''))
bb<-barplot(-log10(DataEachRs_posref$p_wald), names.arg=DataEachRs_posref$rscond, las=3 ,cex.names=0.8, ylab='-log10(pC)')
text(bb[,1], -log10(DataEachRs_posref$p_wald)*0.95,round(DataEachRs_posref$R2,2),cex=0.9)
abline(h=-log10(unique(DataEachRs_posref$p_wald[1])),col='red')
dev.off()

write.csv(DataEachRs_posref, file=paste(opt[['out']],'.csv',sep=''), row.names=F)

