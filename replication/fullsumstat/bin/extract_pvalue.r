#!/usr/bin/env Rscript
library(data.table)
library(qqman)
library("optparse")

option_list = list(
  make_option(c( "--gwas_ref"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_chr"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_bp"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_p"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_a1"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_a2"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_n"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_beta"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_se"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_ref_af"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_chr"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_bp"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_a1"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_a2"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_n"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_af"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_p"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--file_clump"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_beta"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas_cmp_se"), type="character", default=NA,
              help="dataset file name", metavar="character"),
  make_option(c("--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);




getgc<-function(p){
   #https://rdrr.io/cran/gap/src/R/gcontrol2.R
   p <- p[!is.na(p)]
   n <- length(p)
   x2obs <- qchisq(p,1,lower.tail=FALSE)
   x2exp <- qchisq(1:n/n,1,lower.tail=FALSE)
   lambda <- median(x2obs)/median(x2exp)
   return(lambda)
}

getgwas<-function(filegwas, headres,chrgwas, bpgwas, pgwas, a1gwas, a2gwas, betagwas, segwas, n, af){
splline_head<-strsplit(readLines(filegwas,n=1),split='[ \t]')[[1]]
print(splline_head)
data1<-fread(filegwas)
names(data1)<-splline_head
listhead<-c(chrgwas, bpgwas, pgwas, a1gwas, a2gwas, betagwas, segwas,n, af)
listnewhead<-c('chr', 'bp', 'p', 'a1', 'a2', 'beta', 'se', 'n', 'af')
cmthead<-1
for(head in listhead){
 if(!is.na(head)){
  balisehead<-names(data1)==head
  if(!any(balisehead)){
    cat(head,' : ',listnewhead[cmthead] , 'not found')
  }

  cat(head, listnewhead[cmthead],'\n')
  names(data1)[names(data1)==head] <-paste(listnewhead[cmthead], headres,sep='_')
 }
cmthead<-cmthead+1
}
return(data1)
}

#FileSumStat1<-"/home/jeantristan/Travail/GWAS/GWAS_CKD/ImputedDataVAfterReview/GWAS2/All/gemma/all_imputed_map_qc-ckdepi-agesexrank-all.gemma";chrgwas1="chr";bpgwas1="ps";a1gwas1="allele1";a2gwas1="allele0";afgwas1="af";betagwas1="beta";segwas1="se";pgwas1="p_wald";ngwas1=NA;afgwas1='af'
#rsid   chr     bp      a1      a0      beta    se      p       n
#FileSumStat2<-"/home/jeantristan/Data/DataGWASCKD/eGFR/Format/MorrisEtAl2019/COGENT_Kidney_eGFR_trans_ethnic.format.awigen";chrgwas2="chr";bpgwas2="bp";a1gwas2="a1";a2gwas2="a0";afgwas2=NA;betagwas2="beta";segwas2="se";pgwas2="p";ngwas2='n';afgwas2=NA;out<-'teumereGFR_'

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

balisetest<-F
if(balisetest){
opt=list(
gwas_ref='/home/jeantristan/Travail/GWAS/GWAS_CKD/ImputedDataVAfterReview/GWAS2/All/gemma/all_imputed_map_qc-ckdepi-agesexrank-all.gemma', gwas_ref_chr='chr', gwas_ref_bp='ps', gwas_ref_a1='allele1', gwas_ref_a2='allele0', gwas_ref_beta='beta', gwas_ref_se='se', gwas_ref_p='p_wald', gwas_ref_n=NA,gwas_ref_af='af',
#gwas_cmp='/home/jeantristan/Data/DataGWASCKD/eGFR/Format/MorrisEtAl2019/COGENT_Kidney_eGFR_trans_ethnic.format.awigen', gwas_cmp_chr='chr', gwas_cmp_bp='bp', gwas_cmp_a1='a1', gwas_cmp_a2='a0', gwas_cmp_beta='beta', gwas_cmp_se='se', gwas_cmp_p='p', gwas_cmp_n='n',gwas_cmp_af='freqA1',
gwas_cmp='/home/jeantristan/Data/DataGWASCKD/eGFR/Format/MorrisEtAl2019/COGENT_Kidney_eGFR_trans_ethnic.format.awigen', gwas_cmp_chr='chr', gwas_cmp_bp='bp', gwas_cmp_a1='a1', gwas_cmp_a2='a0', gwas_cmp_beta='beta', gwas_cmp_se='se', gwas_cmp_p='p', gwas_cmp_n='n',gwas_cmp_af='freqA1',
out='teumer',
file_clump='ckdepirange_0.000005_250_0.5.clumped'
)
}
filesum_ref<-opt[['gwas_ref']];chr_ref=opt[["gwas_ref_chr"]];bp_ref=opt[["gwas_ref_bp"]];a1_ref=opt[["gwas_ref_a1"]];a2_ref=opt[["gwas_ref_a2"]];beta_ref=opt[['gwas_ref_beta']];se_ref=opt[['gwas_ref_se']];p_ref=opt[['gwas_ref_p']];n_ref=opt[['gwas_ref_n']];af_ref=opt[['gwas_ref_af']]
filesum_cmp<-opt[['gwas_cmp']];chr_cmp=opt[["gwas_cmp_chr"]];bp_cmp=opt[["gwas_cmp_bp"]];a1_cmp=opt[["gwas_cmp_a1"]];a2_cmp=opt[["gwas_cmp_a2"]];beta_cmp=opt[['gwas_cmp_beta']];se_cmp=opt[['gwas_cmp_se']];p_cmp=opt[['gwas_cmp_p']];n_cmp=opt[['gwas_cmp_n']];af_cmp=opt[['gwas_cmp_af']]
headout=opt[['out']]
## clump

#getgwas<-function(filegwas, headres,chrgwas, bpgwas, pgwas, a1gwas, a2gwas, betagwas, segwas, n, af){
dataref<-getgwas(filesum_ref, 'ref', chr_ref, bp_ref, p_ref,a1_ref, a2_ref, beta_ref ,se_ref, n_ref, af_ref)
datacmp<-getgwas(filesum_cmp, 'cmp', chr_cmp, bp_cmp, p_cmp,a1_cmp, a2_cmp, beta_cmp,se_cmp, n_cmp, af_cmp)
dataclump<-read.table(opt[['file_clump']], header=T)

dataclump_ref<-merge(dataclump, dataref, by.x=c('CHR','BP'),by.y=c('chr_ref', 'bp_ref'), all.x=T)
dataclump_all<-merge(dataclump_ref, datacmp, by.x=c('CHR','BP'),by.y=c('chr_cmp', 'bp_cmp'), all.x=T)

balise<-!is.na(dataclump_all$p_cmp)
dataclump_all$p_cmp.adj[balise]<-p.adjust(dataclump_all$p_cmp[balise])
dataclumpsig<-dataclump_all[balise & dataclump_all$p_cmp.adj<0.05,]
write.csv(dataclump_all, file=paste(headout,'_all_clump.csv',sep=''),row.names=F)
write.csv(dataclumpsig, file=paste(headout,'_sig_clump.csv',sep=''),row.names=F)
write.table(dataclumpsig[,c('CHR', 'BP','SNP')], row.names=F, col.names=F, quote=F, file=paste(headout,'_bp',sep=''))



