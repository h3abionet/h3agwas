#!/usr/bin/env Rscript
library(parallel)
library("optparse")
library(data.table)
checkhead<-function(head,Data, type){
if(length(which(head %in% names(Data)))==0){
print(names(Data))
print(paste('not found ', head,'for',type ,'data'))
q(2)
}
}

gopt<-function(x){
gsub('-','.',opt[[x]])
}



option_list = list(
  make_option(c( "--gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--chr_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--ps_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--p_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--wind"), type="numeric", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--nrep"), type="integer", default=NULL,
              help="dataset file name", metavar="integer"),
  make_option(c("--cpus"), type="character", default=NULL,
              help="dataset file name", metavar="integer"),
  make_option(c("--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

headse=gopt('se_gwas');headbp=gopt('ps_gwas');headchr=gopt('chr_gwas');headpval=gopt('p_gwas')
headchrcat=gopt('chr_gwascat');headbpcat=gopt('ps_gwascat');


outhead=opt[['out']]
around=opt[['wind']]*1000
nrep=opt[['nrep']]


datagwas<-fread(opt[['gwas']], header=T)
checkhead(headpval, datagwas,'pval');checkhead(headbp, datagwas,'bp');checkhead(headchr, datagwas, 'chr')


nrep=min(nrep, nrow(datagwas)-1)
listtmp<-sample(1:nrow(datagwas), nrep)
resumeres<-mclapply(listtmp, function(x){
infopos<-as.data.frame(datagwas[x,]);bp<-infopos[,headbp];chr<-infopos[,headchr];begin<-bp-around;end<-bp+around
balise<-(datagwas[[headbp]]>=begin & datagwas[[headbp]]<=end) & datagwas[[headchr]]==chr
pval<-datagwas[balise,..headpval]
Min<-min(pval)
Nb<-nrow(pval)
return(c(Min,Nb))
}, mc.cores=opt[['cpus']])

minpval<-sapply(resumeres, function(x)x[1])
writeLines(as.character(minpval), con=opt[['out']])

