#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--corname"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--data_i"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
   make_option(c("--rel_all"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
   make_option(c("--rel_dup"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
   make_option(c("--missing_dup"), type="character", default=NULL,
              help="output file name [default= %default]", metavar="character"),
   make_option(c("--min_pihat"), type="character", default=NULL,
              help="output file name [default= %default]", metavar="character"),
   make_option(c("--out"), type="character", default=NULL,
              help="head rs [default= %default]", metavar="character")
);

## first steps we identify good related 
args = commandArgs(trailingOnly=TRUE)
if(length(args)<3){
#check_rel.r --data clean_dup.corname  --rel_all out_duplicate.genome --rel_dup out_full.genome  --missing_dup out_duplicate.imiss --min_pihat 0.7 --out out_duplicate
opt<-list(corname='clean_dup.corname', rel_all='out_duplicate.genome', rel_dup='out_full.genome', missing_dup='out_duplicate.imiss', out='missing', min_pihat=0.7, data_i=NULL)
}else{
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
}

dataphenoi<-read.table(opt[['corname']], header=F)
minpihat<-as.numeric(opt[['min_pihat']])

data_reldup<-read.table(opt[['rel_dup']], header=T)
dataphenoi_1<-dataphenoi
dataphenoi_2<-dataphenoi
names(dataphenoi_1)<-c('FID1','IID1', 'FID_N1', 'IID_N1');dataphenoi_1$FID_IID_N1<-paste(dataphenoi_1$FID_N1,dataphenoi_1$IID_N1,sep=':')
names(dataphenoi_2)<-c('FID2','IID2', 'FID_N2', 'IID_N2');dataphenoi_2$FID_IID_N2<-paste(dataphenoi_2$FID_N2,dataphenoi_2$IID_N2,sep=':')
data_reldup<-merge(merge(data_reldup, dataphenoi_1,by=c('FID1','IID1'),all=F), dataphenoi_2, by=c('FID2','IID2'),all=F)
data_reldup2<-data_reldup[data_reldup$FID_IID_N1==data_reldup$FID_IID_N2,]
##bad duplicate 
data_reldup_baddup=data_reldup2[data_reldup2$PI_HAT<minpihat ,]
listbaddup<-as.character(data_reldup_baddup$FID_IID_N1)


datapheno<-dataphenoi
names(datapheno)<-c('FID','IID', 'FID_N', 'IID_N');
datapheno$FID_IID<-paste(datapheno$FID_N, datapheno$IID_N,sep=':')
datapheno$duplicate<-F
tbdup<-table(datapheno$FID_IID)
listdup<-names(tbdup[tbdup>1])
datapheno$duplicate[datapheno$duplicate %in% listdup]<-T

datapheno$bad_duplicate<-F
datapheno$bad_duplicate[datapheno$FID_IID %in% listbaddup]<-T

## inside bad duplicate we researched individual with same pair ()
data_relall<-read.table(opt[['rel_all']], header=T)
data_relall<-merge(merge(data_relall, dataphenoi_1,by=c('FID1','IID1'),all=F), dataphenoi_2, by=c('FID2','IID2'),all=F)
#;data_relall$FID1<-as.character(data_relall$FID1);data_relall$FID2<-as.character(data_relall$FID2)

data_relall_errorother<-data_relall[data_relall$FID_IID_N1!=data_relall$FID_IID_N2 & ((data_relall$FID_IID_N1 %in% listbaddup) | (data_relall$FID_IID_N2 %in% listbaddup)) ,]

#list_badother<-unique(c(data_relall_errorother$FID1, data_relall_errorother$FID2))
#list_badother<-list_badother[!(list_badother %in% listbaddup)]

datapheno$badother_dup<-F
datapheno$badother_dup[datapheno$FID_IID %in% data_relall_errorother$FID_IID_N1]<-T

## check for each 

## we selected individuals with 

## we read missing
datamissing<-read.table(opt[['missing_dup']], header=T)
datamissing<-merge(datamissing,datapheno,by=c('FID','IID'))
datamissing<-datamissing[!(datamissing$FID_IID %in% listbaddup),]
datamissing$FID_IID_I<-paste(datamissing$FID,datamissing$IID,sep=':')
listdup<-unique(as.character(datamissing$FID_IID))
datapheno$FID_IID_I<-paste(datapheno$FID,datapheno$IID,sep=':')
for(dup in listdup){
datamissingdup<-datamissing[datamissing$FID_IID == dup,c('FID','IID','F_MISS', 'FID_IID', 'FID_IID_I')] 
datapheno$dup_choose[datapheno$FID_IID==dup]<-F
lessmiss<-as.character(datamissingdup$FID_IID_I[which.min(datamissingdup$F_MISS)])
datapheno$dup_choose[datapheno$FID_IID_I==lessmiss]<-T
}
#datapheno2<-merge(datapheno,datamissing,all.x=T, by=c('FID','IID'))
datapheno2<-datapheno
##
todelete<-datapheno2$badother_dup | datapheno2$bad_duplicate | ( !is.na(datapheno2$dup_choose) & datapheno2$dup_choose==F)

write.table(datapheno2[todelete, ], file=paste(opt[['out']],'.indtodel.pheno',sep=''), row.names=F, col.names=T, quote=F, sep='\t')
datapheno_clean<-datapheno2[!todelete, ]
write.table(datapheno_clean[, c('FID', 'IID')], file=paste(opt[['out']],'.indtokeep',sep=''), row.names=F, col.names=F, quote=F, sep='\t')

write.table(datapheno_clean[,c('FID','IID', 'FID_N','IID_N')], file=paste(opt[['out']],'.id_update',sep=''),  row.names=F, col.names=F, quote=F,sep='\t')
if(!is.null(opt[['data_i']])){
datapheno_i<-read.table(opt[['data_i']], header=T)
print(head(datapheno_clean))
allpheno<-merge(datapheno_i,datapheno_clean[,c('IID','FID', 'FID_N','IID_N')], by.x=c('FID_update','IID_update'), by.y=c('FID_N','IID_N'),all=F)
allpheno$FID_previous<-allpheno$FID
allpheno$IID_previous<-allpheno$IID
allpheno$FID<-allpheno$FID_update
allpheno$IID<-allpheno$IID_update
tmporder<-names(allpheno)
lishead=c('FID','IID',tmporder[!(tmporder %in% c('FID','IID'))])
write.table(allpheno[,lishead], file=paste(opt[['out']],'.pheno',sep=''),  row.names=F, col.names=T, quote=F,sep='\t')
}

##
#write.table(datapheno_clean, file=paste(opt[['out']],'.full.newpheno',sep=''), col.names=T,row.names=F, quote=F, sep='\t')
#datapheno_clean2<-datapheno_clean[,c('FID', 'IID2', names(datapheno_clean)[!(names(datapheno_clean) %in% c('FID', 'IID', 'IID2'))])]
#names(datapheno_clean2)[2]<-'IID'
#print(table(table(datapheno_clean2[,'FID'])))
#print(grep("_1",datapheno_clean2[,'IID'] ,value=T))
#write.table(datapheno_clean2[,names(dataphenoi)], file=paste(opt[['out']],'.newpheno',sep=''), col.names=T,row.names=F, quote=F, sep='\t')

