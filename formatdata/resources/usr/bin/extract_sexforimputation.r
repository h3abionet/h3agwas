#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## read fam file 
Fam<-read.table(args[1])
InitNAme<-names(Fam)#;Fam$Order<-1:nrow(Fam)
# case of 
headline<-readLines(args[2])
headline<-headline[length(headline)]
Data<-data.frame(Ind=strsplit(headline[length(headline)],split='\t')[[1]]);Data$Ind<-as.character(Data$Ind)


## case 1 : same FID
sameFID<-merge(Data,Fam, by=c(1));nsameFID<-nrow(sameFID)
## case 2 : same ID
sameID<-merge(Data,Fam, by.x=1,by.y=c(2));nsameID<-nrow(sameID)
## case 3 : merge fid and iid
Fam$FIDIID=paste(Fam$V1,Fam$V2,sep='_')
#Data3<-Data;Data3$V1<-sapply(strsplit(as.character(Data3$Ind),split='_'),function(x)x[1]);Data3$V2<-sapply(strsplit(as.character(Data3$Ind),split='_'),function(x)x[2]) 
CompositFIDIID<-merge(Data,Fam, by.x=c('Ind'), by.y='FIDIID');nCompositFIDIID<-nrow(CompositFIDIID)

##case 3 =>  IID : %FID_%IID

if(nsameFID+nsameID+nCompositFIDIID==0){
cat('--------------------------- FAM ----------------')
print(head(Fam[,c(1,2)], 25))
cat('--------------------------- Data ----------------')
print(head(Data[,c(1)], 25))
cat('problem no name corresponding between fam file and data')
quit("no", 5)
}

bestsol<-which.max(c(nsameFID,nsameID,nCompositFIDIID))
print(bestsol)
if(bestsol==1)DataSex<-sameFID
if(bestsol==2)DataSex<-sameID
if(bestsol==3)DataSex<-CompositFIDIID
#5 Sex code ('1' = male, '2' = female, '0' = unknown)
DataSex$Sex<-NA 
DataSex$Sex[DataSex$V5==1]<-'M'
DataSex$Sex[DataSex$V5==2]<-'F'
write.table(DataSex[,c('Ind', 'Sex')], row.names=F, col.names=F, quote=F,sep=' ', file=args[3])
