#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


## read fam file 
Fam<-read.table(args[1], stringsAsFactors=F);Fam<-Fam[,c(1,2)];names(Fam)<-c('FID','IID');Fam$FID_I<-Fam$FID;Fam$IID_I<-Fam$IID
InitNAme<-names(Fam)#;Fam$Order<-1:nrow(Fam)
Data<-read.table(args[2], header=T, stringsAsFactors=F)[,c(1,2)]

## case 1 
sameFIDIID<-merge(Data,Fam, by=c(1,2), all=F);nsameFIDIID<-nrow(sameFIDIID)

##case 2 : information just on ID  (FID)
Fam2<-Fam;Fam2$FID<-Fam2$IID
sameIID<-merge(Data,Fam2, by=c(1,2), all=F);nsameIID<-nrow(sameIID)
head(sameIID)

##case 3 =>  IID : %FID_%IID
Fam3<-Fam;Fam3$FID<-sapply(strsplit(as.character(Fam2$IID),split='_'),function(x)x[1]);Fam3$IID<-sapply(strsplit(as.character(Fam2$IID),split='_'),function(x)x[2]) 
CompositFIDIID<-merge(Data,Fam3, by=c(1,2), all=F);nCompositFIDIID<-nrow(CompositFIDIID)
head(CompositFIDIID)

if(nsameFIDIID+nsameIID+nCompositFIDIID==0){
cat('--------------------------- FAM ----------------')
print(head(Fam[,c(1,2)], 25))
cat('--------------------------- Data ----------------')
print(head(Data[,c(1,2)], 25))
cat('problem no name corresponding between fam file and data')
quit("no", 5)
}

bestsol<-which.max(c(nsameFIDIID,nsameIID,nCompositFIDIID))
head(CompositFIDIID)
if(bestsol==1)write.table(sameFIDIID[,c('FID_I','IID_I', 'FID','IID')], file=args[3], row.names=F, col.names=F, quote=F, sep='\t')
if(bestsol==2)write.table(SameIID[,c('FID_I','IID_I', 'FID','IID')], file=args[3], row.names=F, col.names=F, quote=F, sep='\t')
if(bestsol==3)write.table(CompositFIDIID[,c('FID_I','IID_I', 'FID','IID')], file=args[3], row.names=F, col.names=F, quote=F, sep='\t')
