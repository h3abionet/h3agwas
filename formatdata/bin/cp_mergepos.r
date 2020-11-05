#!/usr/bin/env Rscript
library("optparse")
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
   make_option(c("-z", "--file_rsres"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
   make_option(c("-w", "--file_cross"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
   make_option(c("-o", "--out"), type="character", default=NULL,
              help="output file name [default= %default]", metavar="character"),
   make_option(c("-r", "--head_rs"), type="character", default=NULL,
              help="head rs [default= %default]", metavar="character"),
   make_option(c("-b", "--head_bp"), type="character", default=NULL,
              help="head rs [default= %default]", metavar="character"),
   make_option(c("-c", "--head_chr"), type="character", default=NULL,
              help="head rs [default= %default]", metavar="character"),
   make_option(c("-s", "--sep"), type="character", default=NULL,
              help="head rs [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser);


FileI=args[['file']]
FileRsRes=args[['file_rsres']]
FileRsCros=args[['file_cross']]
FileOutRs=paste(args[['out']],'.rs',sep='')
FileOutPos=paste(args[['out']],'.pos',sep='')
Sep=GetSep(args[['sep']])
ChrHead=args[['head_chr']]
BpHead=args[['head_bp']]
RsHead=args[['head_rs']]


Data<-read.table(FileI, header=T,sep=Sep,comment.char="", quote="", stringsAsFactors=F)
HeaderI<-strsplit(gsub('\n','',readLines(FileI, 1)),split=Sep)[[1]]
names(Data)<-HeaderI

PosBegin<-sapply(strsplit(as.character(Data[,BpHead]),split=";"),function(x){
min(as.integer(x,na.rm=T))
})
PosEnd<-sapply(strsplit(as.character(Data[,BpHead]),split=";"),function(x){
max(as.integer(x), na.rm=T)
})
Chr<-sapply(strsplit(as.character(Data[,ChrHead]),split=";"),function(x){
Un<-unique(x[as.integer(x)>0])
if(length(Un)==1)return(Un)
else return(NA)
}
)

Data$PosBeginI<-PosBegin
Data$PosEndI<-PosEnd
Data$ChrI<-Chr
Data$Num<-1:nrow(Data)

LineCrossMap<-readLines(FileRsCros);nbcol<-sapply(strsplit(LineCrossMap, split="\t"),function(x)length(x))
LineCrossMap<-LineCrossMap[nbcol==5]
writeLines(LineCrossMap, con='tmplaflkjv')
DataCrossMap<-read.table('tmplaflkjv',header=F)
names(DataCrossMap)<-c("ChroNewCM","PosDebNewCM","PosFinNewCM",'Strand' ,"Num")

DataRs<-read.table(FileRsRes,sep="\t")
DataRs<-DataRs[,c("V1","V2","V3","V4","V5")]
names(DataRs)<-c("ChroNewRs","PosNewRs","Rs","Ref","Alt")

DataM<-merge(merge(Data,DataRs,by.x="SNPS",by.y="Rs", all=T),DataCrossMap, by.x="Num", by.y="Num", all=T)
DataNotFound<-DataM[(is.na(DataM$ChroNewCM) & is.na(DataM$ChroNewRs)),]
write.table(DataNotFound, file=paste(args[['out']],".notfound.tsv", sep='') ,sep='\t', row.names=F, col.names=T)
DataM<-DataM[!(is.na(DataM$ChroNewCM) & is.na(DataM$ChroNewRs)),]
write.table(DataM, file=paste(args[['out']],".detail.tsv", sep='') ,row.names=F, sep='\t',  col.names=T)

DataM$ChroNew<-as.character(DataM$ChroNewRs)
DataM$PosBeginNew<-DataM$PosNewRs
DataM$PosEndNew<-DataM$PosNewRs
Bal<-is.na(DataM$ChroNew)
DataM$ChroNew[Bal]<-DataM$ChroNewCM[Bal]
DataM$PosBeginNew[Bal]<-DataM$PosDebNewCM[Bal]
DataM$PosEndNew[Bal]<-DataM$PosFinNewCM[Bal]

DataM$ChroNew<-gsub("chr","",DataM$ChroNew)
tmpsup1<-table(DataM$Num)
tmpsup1<-as.integer(names(tmpsup1[tmpsup1>1]))
DataMMulti<-DataM[DataM$Num %in% tmpsup1,]
DataM<-DataM[!(DataM$Num %in% tmpsup1),]
DataM<-DataM[!is.na(DataM$ChroNew) ,!(names(DataM) %in% c("Num", "ChroNewCM","PosDebNewCM","PosFinNewCM","ChroNewRs","PosNewRs"))]

write.table(DataM , sep="\t", row.names=F,col.names=T,file=paste(args[['out']],".tsv", sep='') )
write.table(DataMMulti, sep="\t", row.names=F,col.names=T,file=paste(args[['out']],".multi.tsv", sep='') )


