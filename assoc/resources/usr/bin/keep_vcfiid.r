#!/usr/bin/env Rscript
library("optparse")
#     addpheno_pcs.r --data $phenofile --pcs $head".eigenvec" --out newfilepheno
option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--vcf"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt= parse_args(opt_parser);


## read fam file 

Data<-read.table(opt[['data']], stringsAsFactors=F);print(head(Data));Data<-Data[,c(1,2)];names(Data)<-c('FID','IID');Data$FID_I<-Data$FID;Data$IID_I<-Data$IID

vcf_name<-readLines(opt[['vcf']])
## case 1 
baliseFID<-vcf_name %in% as.character(Data$FID);nFID<-length(which(baliseFID))
baliseIID<-vcf_name %in% as.character(Data$IID);nIID<-length(which(baliseIID))
FIDIID<-paste(as.character(Data$FID), as.character(Data$IID),sep='_')
baliseFIDIID<-vcf_name %in% FIDIID;nFIDIID<-length(which(baliseFIDIID))
##case 2 : information just on ID  (FID)

##case 3 =>  IID : %FID_%IID
if(nIID+nFID+nFIDIID==0){
cat('--------------------------- FAM ----------------')
print(head(Data[,c(1,2)], 25))
cat('--------------------------- Data ----------------')
print(head(vcfname, 25))
cat('problem no name corresponding between vcf and data')
quit("no", 5)
}

bestsol<-which.max(c(nIID,nFID,nFIDIID))
if(bestsol==1)writeLines(vcf_name[baliseIID], con=opt[['out']])
if(bestsol==2)writeLines(vcf_name[baliseFID], con=opt[['out']])
if(bestsol==3)writeLines(vcf_name[baliseFIDIID], con=opt[['out']])
