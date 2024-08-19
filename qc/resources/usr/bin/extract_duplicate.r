#!/usr/bin/env Rscript
library("optparse")
checkhead<-function(x, headfile, file){
if(any(!(x %in% headfile))){
cat(x[!(x %in% headfile)], ' not found in ',file,'\nexit\n')
exit(3)
}
}

#   extract_duplicate.r --data $data --bfile $bfile --out $out --col_fidid ${params.col_fidid}  --col_newfidid ${params.newcol_fidid}
option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--bfile"), type="character",
              help="dataset file name", metavar="character"),
  make_option(c("--col_fidid"), type="character", default="FID,IID",
              help="dataset file name", metavar="character"),
  make_option(c("--col_newfidid"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--out"), type="character", default=NULL,
              help="dataset file name", metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
out=opt[['out']]
#opt=list(data="allcancer.changename.tsv",bfile="allind_tmp",out="clean_dup",col_fidid="FID,IID",col_newfidid="IID_update")

Data<-unique(read.table(opt[['data']],header=T))
fam<-read.table(paste(opt[['bfile']], '.fam',sep=''))
fam$V1<-as.character(fam$V1)
fam$V2<-as.character(fam$V2)

splfid=strsplit(opt[['col_fidid']],split=',')[[1]]
checkhead(splfid, names(Data), opt[['data']])
splnewfid=strsplit(opt[['col_newfidid']],split=',')[[1]]
checkhead(splnewfid, names(Data), opt[['data']])
NewData<-Data

NewData$FID<-Data[,splfid[1]]
if(length(splfid)==1)NewData$IID<-Data[,splfid[1]] else {
NewData$IID<-Data[,splfid[2]]
}

NewData$FID_update<-Data[,splnewfid[1]]
if(length(splnewfid)==1)NewData$IID_update<-Data[,splnewfid[1]] else {
NewData$IID_update<-Data[,splnewfid[2]]
}
good<-paste(NewData$FID,NewData$IID) %in% paste(fam$V1,fam$V2)
if(length(which(good))<2){
print(paste(fam$V1,fam$V2))
print(paste(NewData$FID,NewData$IID))
cat('error between bfile and phenotype file\nexit\n')
exit(4)
}
error<-NewData[!good,]
write.csv(error, file=paste(out,'_phenoind_notfound.csv',sep=''),row.names=F)

good_fam<- paste(fam$V1,fam$V2) %in% paste(NewData$FID,NewData$IID) 
errorfam<-fam[!good_fam,]
write.csv(errorfam, file=paste(out,'_fam_notfound.csv',sep=''), row.names=F)

NewData<-NewData[good,]

iidupdate<-paste(NewData$FID_update,NewData$IID_update)
tbiidupdate<-table(iidupdate)
iddup<-names(tbiidupdate[tbiidupdate>1])
NewData_dup<-NewData[iidupdate %in% iddup,c('FID','IID', 'FID_update','IID_update')]
write.table(NewData[,c('FID','IID', 'FID_update','IID_update')], file=paste(out,'.corname',sep=''), row.names=F, col.names=F, sep='\t', quote=F)
write.table(NewData, file=paste(out,'.update',sep=''), row.names=F, col.names=T, sep='\t', quote=F)
if(nrow(NewData[iidupdate %in% iddup,])==0){
cat('no individuals duplicate exit\n')
quit(save = "no", status = 1, runLast = FALSE)
}
write.table(NewData[iidupdate %in% iddup,c('FID','IID', 'FID_update','IID_update')], file=paste(out,'.dup',sep=''), row.names=F, col.names=F, sep='\t',quote=F)
