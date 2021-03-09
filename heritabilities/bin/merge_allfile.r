#!/usr/bin/env Rscript
library(ggplot2)
args = commandArgs(trailingOnly=TRUE)


if (length(args)!=2){
  stop("At least two argument must be supplied input file and models", call.=FALSE)
}

listfile=strsplit(args[1], split=',')[[1]]
listfile="1-pheno-1_type1_gemma.stat"
headout=args[2]
Cmt<-1
for(file in listfile){
Data<-read.table(file, header=T)
if(Cmt==1)DataF<-Data
else DataF<-rbind(DataF, Data)
Cmt<-Cmt+1
}
write.csv(DataF, row.names=F, file=paste(headout,'.csv',sep=''))

DataF$InfoSoft[DataF$InfoSoft=="None"]<-'-'
DataF$Info<-gsub('^[0-9]+-','', DataF$Info)
DataF$Header<-paste( DataF$Info,',',DataF$TypeCor,' ',DataF$InfoSoft,'',sep='')

#mWg0<-max(resum2NoG0$Perc)
DataF$Header<-as.factor(DataF$Header)
#ggplot(data=DataF, aes(x=Gen, y=reorder(Header, order))) +
ggplot(data=DataF, aes(x=Header, y=Gen)) +
geom_bar(stat="identity", position=position_dodge())+
  labs(x="", y = "Heritabilities") +scale_fill_brewer(palette="Blues") +  coord_flip() #+   scale_fill_manual(values=c('#00C0B8', '#00A5FF', '#F8766D'))
ggsave(paste(headout,'.svg',sep=''), width = 7, height = 7*2)

