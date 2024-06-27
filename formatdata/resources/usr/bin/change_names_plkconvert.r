#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


## read fam file 
Fam<-read.table(args[1])
InitNAme<-names(Fam)#;Fam$Order<-1:nrow(Fam)
Data<-read.table(args[2], header=T)[,c(1,2)]

## case 1 
sameFIDIID<-merge(Data,Fam, by=c(1,2));nsameFIDIID<-nrow(sameFIDIID)

##case 2 : information just on ID  (FID)
Fam2<-Fam;Fam2$V1<-Fam2$V2
sameIID<-merge(Data,Fam2, by=c(1,2));nsameIID<-nrow(sameIID)

##case 3 =>  IID : %FID_%IID
Fam3<-Fam;Fam3$V1<-sapply(strsplit(as.character(Fam2$V2),split='_'),function(x)x[1]);Fam3$V2<-sapply(strsplit(as.character(Fam2$V2),split='_'),function(x)x[2]) 
CompositFIDIID<-merge(Data,Fam3, by=c(1,2));nCompositFIDIID<-nrow(CompositFIDIID)

if(nsameFIDIID+nsameIID+nCompositFIDIID==0){
cat('--------------------------- FAM ----------------')
print(head(Fam[,c(1,2)], 25))
cat('--------------------------- Data ----------------')
print(head(Data[,c(1,2)], 25))
cat('problem no name corresponding between fam file and data')
quit("no", 5)
}

bestsol<-which.max(c(nsameFIDIID,nsameIID,nCompositFIDIID))
if(bestsol==2)write.table(Fam2, file=args[3], row.names=F, col.names=F, quote=F, sep='\t')
if(bestsol==3)write.table(Fam3, file=args[3], row.names=F, col.names=F, quote=F, sep='\t')
