#!/usr/bin/env Rscript
library(data.table)
library("optparse")
library("qqman")
t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color

## Get RGB values for named color
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
             max = 255,
             alpha = (100 - percent) * 255 / 100,
             names = name)

## Save the color
invisible(t.col)
}



computedher<-function(beta, se, af,N){
#https://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0120758.s001
maf<-af
maf[!is.na(af) & af>0.5]<- 1 - maf[!is.na(af) & af>0.5]
ba<-!is.na(beta) & !is.na(se) & !is.na(maf) & !is.na(N)
a<-rep(NA, length(beta))
b<-rep(NA, length(beta))
a<-2*(beta[ba]**2)*(maf[ba]*(1-maf[ba]))
b<-2*(se[ba]**2)*N[ba]*maf[ba]*(1-maf[ba])
res<-rep(NA, length(beta))
res[ba]<-a[ba]/(a[ba]+b[ba])
res
}



option_list = list(
  make_option(c( "--gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c( "--gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c( "--clump_file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c( "--pheno"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--chr_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--bp_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--ps_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--a1_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--a2_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--bim"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--beta_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--se_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--af_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--N_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--ps_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--chr_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--p_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--min_pvalue"), type="numeric", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--min_r2"), type="numeric", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--info_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

checkhead<-function(head,Data, type){
if(length(which(head %in% names(Data)))==0){
print(names(Data))
print(paste('not found ', head,'for',type ,'data'))
q(2)
}
}

#--chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#opt=list(gwascat='out_all.csv',gwas='out_range.init',chr_gwas='chr',ps_gwas='ps',a1_gwas='allele1',a2_gwas='allele0',beta_gwas='beta',se_gwas='se',af_gwas='af',chr_gwascat='chrom',bp_gwascat='chromEnd',p_gwas='p_wald',ps_gwascat='chromEnd',chr_gwascat='chrom','out'='out_ld',clump_file='out.clumped',min_pvalue=0.1,min_r2=0.1,info_gwascat="pubMedID;author;trait;initSample",bim='all_imputed_map_qc.bim')

headse=opt[['se_gwas']];headbp=opt[['ps_gwas']];headchr=opt[['chr_gwas']];headbeta=opt[['beta_gwas']];heada1=opt[['a1_gwas']];heada2=opt[['a2_gwas']];headpval=opt[['p_gwas']];headaf<-opt[['af_gwas']];headbeta=opt[['beta_gwas']]
headchrcat=opt[['chr_gwascat']];headbpcat=opt[['ps_gwascat']];heada1catrs<-"riskAllele";headzcat="z.cat";headafcat<-'risk.allele.af';heada1cat<-'risk.allele.cat'
outhead=opt[['out']]

filegwascat<-opt[['gwascat']];filegwas=opt[['gwas']];fileclump=opt[['clump_file']];filebim=opt[['bim']]
#filegwascat<-"out_all.csv";filegwas="out_pos.init";fileclump="out.clumped";filebim="all_imputed_map_qc.bim"



datagwascat=read.csv(filegwascat)
datagwascat[,heada1cat]<-sapply(strsplit(as.character(datagwascat[,heada1catrs]),split='-'),function(x)x[2])
datagwas<-read.table(filegwas, header=T)
checkhead(headaf, datagwas,'af');checkhead(headpval, datagwas,'pval');checkhead(headse, datagwas,'se');checkhead(headbp, datagwas,'bp');checkhead(headchr, datagwas, 'chr');checkhead(headbeta, datagwas, 'beta')

checkhead(headbpcat,datagwascat,'bp cat');checkhead(headchrcat,datagwascat,'chro cat');

# CHR_A         BP_A         SNP_A  CHR_B         BP_B         SNP_B           R2 
dataclump<-read.table(fileclump,header=T)
databim<-fread(filebim)[,c(1,4,2)];names(databim)<-c('chr','bp','rs_bim')

datagwascat<-merge(datagwascat, databim,by.x=c(headchrcat, headbpcat), by.y=c('chr','bp'),all=F)
datagwas<-merge(datagwas, databim,by.x=c(headchr, headbp), by.y=c('chr','bp'),all=F)





headN<-opt[['N_value']]
if(is.null(opt[['N_gwas']])){
if(is.null(opt[['N_value']]))Nval<-10000 else Nval=opt[['N_value']]
datagwas[,'N_gwas']<-Nval
headN<-'N_gwas'
}


datagwas$h2.gwas<-computedher(datagwas[,headbeta], datagwas[,headse], datagwas[,headaf],datagwas[,headN])
datagwas$z.gwas<-datagwas[,headbeta]/datagwas[,headse]
datagwasf<-datagwas[datagwas[,headpval]<headpval,]

tmpsp2<-strsplit(gsub('\\([1-9]+\\)','',as.character(dataclump$SP2)), split=',')
names(tmpsp2)<-as.character(dataclump$SNP)
cmt<-1
for(wind in names(tmpsp2)){
tmprs<-data.frame(rsclump=wind, rs_wind=c(tmpsp2[[wind]],wind))
if(cmt==1)resclump<-tmprs
else resclump<-rbind(resclump ,tmprs)
cmt<-cmt+1
}
resclumpgwas<-merge(resclump,datagwas,by.x='rs_wind',by.y="rs_bim")
resclumpgwascat<-merge(resclump,datagwascat,by.x='rs_wind',by.y="rs_bim")
resall<-merge(resclumpgwas,resclumpgwascat,by=c('rsclump'), all=T,suffixes = c("_gwas","_gwascat"))
write.csv(resall,file=paste(opt[['out']],'_detailall.csv',sep=''), row.names=F)
balise<-!is.na(resall$rs_wind_gwas)
resall$info_gwas[balise]<-paste(resall[balise,headchr],':',resall[balise,headbp],'-beta:',resall[balise,headbeta], ',se:',resall[balise,headse],',pval:',resall[balise,headpval],sep='')

resall$info_gwascat<-""
infocat=strsplit(opt[['info_gwascat']],split=';')[[1]]
balise<-!is.na(resall$rs_wind_gwascat)
if(all(balise==F)){
datacatinfo<-resall[F,c('info_gwascat','rsclump')]
datacatinfors<-resall[F,c('rs_wind_gwascat','rsclump')]
datainfo<-resall[F,c('info_gwas','rsclump')]
datainfors<-resall[F,c('rs_wind_gwas','rsclump')]
}else{
resall$info_gwascat[balise]=paste(resall$rs_wind_gwascat[balise],':',resall[balise,headchrcat], ':',resall[balise,headbpcat],',',sep='')
for(cat in infocat)resall$info_gwascat[balise]<-paste(resall$info_gwascat[balise],cat,':',resall[balise,cat],',',sep='')

datacatinfo<-aggregate(info_gwascat~rsclump,data=resall,FUN=function(x)paste(unique(x),collapse=';'))
datacatinfors<-aggregate(rs_wind_gwascat~rsclump,data=resall,FUN=function(x)paste(unique(x),collapse=';'))

datainfo<-aggregate(info_gwas~rsclump,data=resall,FUN=function(x)paste(unique(x),collapse=';'))
datainfors<-aggregate(rs_wind_gwas~rsclump,data=resall,FUN=function(x)paste(unique(x),collapse=';'))
}

allcati<-merge(merge(merge(datacatinfo,datacatinfors, by='rsclump'), datainfo,by='rsclump'),datainfors,by='rsclump')
allcati<-merge(allcati,dataclump[,c('SNP','CHR','BP','TOTAL','NSIG', 'P')],by.x='rsclump', by.y='SNP',all=F)
write.csv(allcati ,file=paste(opt[['out']],'_allresum.csv',sep=''), row.names=F)

