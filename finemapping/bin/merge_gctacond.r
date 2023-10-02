#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("--files_cond"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--ld"), type="character", default="out.txt",
              help="out prefix", metavar="character"),
  make_option(c("--pos_ref"), type="integer", 
              help="out prefix", metavar="character"),
  make_option(c("--out"), type="character",
              help="list of genes ", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

pos_ref=as.integer(opt[['pos_ref']])

list_cond=strsplit(opt[['files_cond']],split=',')[[1]]
Cmt<-1
for(File in list_cond){
 rsinfo=strsplit(File,split='_')[[1]][1]
 chro=strsplit(rsinfo,split=':')[[1]]
 bp=chro[1]
 Data<-read.table(File,header=T)
 if(!is.null(pos_ref)){
  Data<-Data[Data$bp==pos_ref,]
 }
 Data$refA[as.logical(Data$refA)]<-"T"
 Data$bp_cond<-bp
 Data$bp_cond<-as.character(Data$bp_cond)
 if(Cmt==1)DataF<-Data
 else DataF<-rbind(DataF,Data)
 Cmt<-Cmt+1
}

DataLD<-read.table(opt[['ld']],header=T)
DataLD2<-DataLD
DataLD2$CHR_A<-DataLD$CHR_B;DataLD2$BP_A<-DataLD$BP_B;DataLD2$SNP_A<-DataLD$SNP_B
DataLD2$CHR_B<-DataLD$CHR_A;DataLD2$BP_B<-DataLD$BP_A;DataLD2$SNP_B<-DataLD$SNP_A
DataLD<-rbind(DataLD,DataLD2)
names(DataF)[2]<-'SNP_ref'
names(DataF)[3]<-'bp_ref'

DataF$bp_ref<-as.character(DataF$bp_ref)
DataF$bp_cond<-as.character(DataF$bp_cond)

DataLD$BP_A<-as.character(DataLD$BP_A);DataLD$BP_B<-as.character(DataLD$BP_B);

DataAll<-merge(DataF,DataLD[,c('BP_A','BP_B','R2')], by.x=c('bp_ref','bp_cond'),by.y=c('BP_A','BP_B'),all.x=T)

DataAll<-DataAll[,c("Chr",'bp_cond',"bp_ref","SNP_ref","freq","b","se","p","n","freq_geno","bC","bC_se","pC","R2")]
newnames<-c("Chr",'bp_cond',"bp_ref","SNP_ref","freq_ref","b_ref","se_ref","p_ref","n_ref","freq_geno_ref","bC_ref","bC_se_ref","pC_ref","R2")
names(DataAll)<-newnames

write.csv(DataAll, file=opt[['out']], row.names=F)


