#!/usr/bin/env Rscript
library("optparse")
library(data.table)
library(fastman)

option_list = list(
  make_option(c( "--data"), type="character",
              help="data summary statistics", metavar="character"),
  make_option(c("--p_header"), type="character",
              help="ped file contains genotype for each individual", metavar="character"),
  make_option(c( "--chr_header"), type="character",
              help="phenotype in data", metavar="character"),
  make_option(c( "--bp_header"), type="character",
              help="phenotype in data", metavar="character"),
  make_option(c( "--af_header"), type="character",
              help="phenotype in data", metavar="character"),
  make_option(c( "--rs_header"), type="character",
              help="phenotype in data", metavar="character"),
  make_option(c( "--maf"), 
              help="phenotype in data", type="numeric", default=0.01),
  make_option(c( "--file_rs"), 
              help="phenotype in data", type="character"),
  make_option(c( "--type_plot"), 
              help="phenotype in data", type="character", default='png'),
  make_option(c( "--out"), type="character",
              help="phenotype in data", metavar="character")
);

opt = OptionParser(option_list=option_list);
opt = parse_args(opt);
test=F
if(test){
opt=list("data"="gwas_merge/dbp/env12/dbp_all_env12.gwas","af_header"="EAF_ALL", "p_header"="P_INT_ROBUST","chr_header"="CHR",out= "gwas_merge/dbp/env12/dbp_all_env12")
}

alldata<-fread(opt[['data']])
listhead<-c(opt[['chr_header']], opt[['bp_header']], opt[['p_header']], opt[['rs_header']], opt[['af_header']])
if(all(listhead %in% names(alldata))==F){
cat('header not found', paste(listhead[!(listhead %in% names(alldata))], collapse=','),'\nexit\n')
q('n')
}

out=opt[['out']]
type_plot=opt[['type_plot']]
plotf= get(type_plot)
width=   480
height=   480

if(type_plot  %in% c('pdf','svg') ) width=7 
if(type_plot %in% c('pdf','svg')) height=7

plotf(paste(out,'_distchro.',type_plot,sep=''),   width = width*3, height = height*2)
plot(table(alldata[[opt[['chr_header']]]]), xlab='chromosome', ylab='distribution')
dev.off()
if(!is.null(opt[['af_header']])){
 maf=opt[['maf']];maf2<-1 - maf
 alldataqq<-alldata[alldata[[opt[['af_header']]]]>maf & alldata[[opt[['af_header']]]]<maf2,]
}else {
alldataqq<-alldata
}
alldataqq<-na.omit(alldataqq[,..listhead])

rsheader=opt[['rs_header']]
if(!is.null(opt[['file_rs']])){
data_rs<-fread(opt[['file_rs']],header=F);names(data_rs)<-c('CHR','BP','A1','A2','SNP','REF')
if('SNP' %in% names(alldataqq)){
 alldataqq <- alldataqq[, -which(names(alldataqq) == "SNP" )]
}
alldataqq<-merge(alldataqq,data_rs,by.y=c('CHR','BP'), by.x=c(opt[['chr_header']], opt[['bp_header']]),all.x=T)
rsheader='SNP'
}

print(alldataqq)

plotf(paste(out,'_qq.',type_plot,sep=''),width = width, height = height)
fastqq (alldataqq[[opt[['p_header']]]])
dev.off()

if(is.null(rsheader)){
alldataqq$rsid<-paste(alldataqq[[opt[['chr_header']]]], alldataqq[[opt[['bp_header']]]], sep=':')
rsheader='rsid'
}else{
balise<-is.na(alldataqq[[rsheader]]) | alldataqq[[rsheader]]=='.' 
alldataqq[[rsheader]][balise]<-paste(alldataqq[[opt[['chr_header']]]][balise], alldataqq[[opt[['bp_header']]]][balise], sep=':')
}
if(type_plot  %in% c('pdf','svg'))plotf(paste(out,'_man.',type_plot, sep=''),   width = width*6, height = height*2) else plotf(paste(out,'_man.',type_plot, sep=''),   width = width*6, height = height*2, res=200)
fastman(alldataqq,  chr =opt[['chr_header']] , bp = opt[['bp_header']], p = opt[['p_header']], snp=rsheader,annotatePval=5E-8)
dev.off()


