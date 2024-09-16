#!/usr/bin/env Rscript
library("optparse")
library(data.table)
GetSep<-function(x){
listsep=c(',', ' ', '\t')
listsepnum<-c('COM', 'SPA', 'TAB')
if(x %in% listsep)return(x)
x<-toupper(x)
x<-substr(x,1,3)
if(x %in% listsepnum)return(listsep[listsepnum==x])
cat('\nnot found sep ', x, '\nexit\n')
q(save='no',status=1)
}

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
   make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character"),
   make_option(c("-r", "--head_rs"), type="character", default="out.txt", 
              help="head rs [default= %default]", metavar="character"),
   make_option(c("-b", "--head_bp"), type="character", default="out.txt", 
              help="head rs [default= %default]", metavar="character"),
   make_option(c("-c", "--head_chr"), type="character", default="out.txt", 
              help="head rs [default= %default]", metavar="character"),
   make_option(c("-s", "--sep"), type="character", default="out.txt", 
              help="head rs [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);


FileI=args[['file']]
FileOutRs=paste(args[['out']],'.rs',sep='')
FileOutPos=paste(args[['out']],'.pos',sep='')
Sep=GetSep(args[['sep']])
cat(Sep)
ChrHead=args[['head_chr']]
BpHead=args[['head_bp']]
RsHead=args[['head_rs']]

#Data<-read.table('gwas_catalog.tsv', header=T, sep='\t', comment.char="",quote="")
#Data<-read.table(FileI, header=T,sep=Sep,comment.char="", quote="", stringsAsFactors=F)
Data<-fread(FileI, header=T,sep=Sep, quote="", stringsAsFactors=F)
print(head(Data))
headall<-c(ChrHead,BpHead,RsHead)
Data2<-Data[, ..headall];names(Data2)<-c('Chr', 'Pos', 'rsid')
print(names(Data2))
Data$Chr<-Data[[ChrHead]]
Data$PosBegin<-Data[[BpHead]]
Data$PosEnd<-Data[[BpHead]]
Data$rsid<-Data[[RsHead]]
tmp2<-c('Chr', 'PosBegin', 'PosEnd', 'rsid')
Tmp<-Data[ , ..tmp2]
ListRs<-unlist(unlist(strsplit(Tmp[['rsid']],split=';')))
writeLines(ListRs, con=FileOutRs)
#Tmp$Chr<-as.character(Tmp$Chr)
#bal=grep('chr',Tmp$Chr, invert=T)
#Tmp$Chr[bal]<-paste("chr",as.character(Tmp$Chr[bal]),sep="")

#write.table(unique(Tmp), quote=F, row.names=F, col.names=F,file=FileOut)

#writeLines(unlist(as.vector(strsplit(as.character(Data$SNPS),split="[,;]"))), con="rsTosearch")
PosBegin<-sapply(strsplit(as.character(Data2$Pos),split=";"),function(x){
min(as.integer(x,na.rm=T))
})
PosEnd<-sapply(strsplit(as.character(Data2$Pos),split=";"),function(x){
max(as.integer(x), na.rm=T)
})
Chr<-sapply(strsplit(as.character(Data2$Chr),split=";"),function(x){
Un<-unique(x[as.integer(gsub('chr','',x))>0])
if(length(Un)==1)return(Un)
else return(NA)
}
)
Data2$PosBeginI<-PosBegin
Data2$PosEndI<-PosEnd
Data2$ChrI<-Chr
Data2$RsId<-Data2$rsid
Data2$Num<-1:nrow(Data2)
Bal<- !is.na(Data2$Chr) & !is.na(Data2$PosBegin) & !is.infinite(Data2$PosBegin)
bal=grep('chr',Data2$ChrI, invert=T)
Data2$ChrI[bal]<-paste("chr",as.character(Data2$ChrI[bal]),sep="")
Data2$Strand<-'+'
Data2ToP<-Data2[Bal, c("ChrI","PosBeginI","PosEndI", 'Strand',"Num")]

write.table(Data2ToP, file=FileOutPos, sep="\t", quote=F, row.names=F, col.names=F)



