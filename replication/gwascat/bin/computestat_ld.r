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
  make_option(c( "--ld_file"), type="character", default=NULL,
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
Test=F
#--chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if(Test)opt=list(gwascat='Diastolic_AwigenLD_all.csv',gwas='Diastolic_AwigenLD_range.init',chr_gwas='CHR',ps_gwas='BP',a1_gwas='ALLELE1',a2_gwas='ALLELE0',beta_gwas='BETA',se_gwas='SE',af_gwas='A1FREQ',chr_gwascat='chrom',bp_gwascat='chromEnd',p_gwas='P_BOLT_LMM',ps_gwascat='chromEnd',chr_gwascat='chrom',out='Diastolic_AwigenLD_ld',ld_file='Diastolic_AwigenLD_ld.ld',min_pvalue='0.001',min_r2=0.2,info_gwascat="pubMedID;author;trait;initSample")


headse=opt[['se_gwas']];headbp=opt[['ps_gwas']];headchr=opt[['chr_gwas']];headbeta=opt[['beta_gwas']];heada1=opt[['a1_gwas']];heada2=opt[['a2_gwas']];headpval=opt[['p_gwas']];headaf<-opt[['af_gwas']];headbeta=opt[['beta_gwas']]
headchrcat=opt[['chr_gwascat']];headbpcat=opt[['ps_gwascat']];heada1catrs<-"riskAllele";headzcat="z.cat";headafcat<-'risk.allele.af';heada1cat<-'risk.allele.cat'
outhead=opt[['out']]


datagwascat=read.csv(opt[['gwascat']])
#datagwascat[,heada1cat]<-sapply(strsplit(as.character(datagwascat[,heada1catrs]),split='-'),function(x)x[2])
datagwas<-read.table(opt[['gwas']], header=T)
checkhead(headaf, datagwas,'af');checkhead(headpval, datagwas,'pval');checkhead(headse, datagwas,'se');checkhead(headbp, datagwas,'bp');checkhead(headchr, datagwas, 'chr');checkhead(headbeta, datagwas, 'beta')

checkhead(headbpcat,datagwascat,'bp cat');checkhead(headchrcat,datagwascat,'chro cat');

# CHR_A         BP_A         SNP_A  CHR_B         BP_B         SNP_B           R2 
datald<-fread(opt[['ld_file']])
tmpa<-unique(datald[,c('CHR_A','BP_A','SNP_A')]);names(tmpa)<-c('CHR','BP','SNP');
tmpb<-unique(datald[,c('CHR_B','BP_B','SNP_B')]);names(tmpb)<-c('CHR','BP','SNP');
tmp<-rbind(tmpa,tmpb)
tmpall<-unique(rbind(tmpa,tmpb))
tmpall<-tmpall[,c('CHR','BP','SNP','CHR','BP','SNP')]
tmpall$R2<-1
names(tmpall)<-names(datald)
datald<-rbind(datald,tmpall)

#   CHR_A     BP_A     SNP_A CHR_B     BP_B      SNP_B       R2
#1:    18 48132646 rs1437649    18 48133241 rs61148001 0.896039





headN<-opt[['N_value']]
if(is.null(opt[['N_gwas']])){
if(is.null(opt[['N_value']]))Nval<-10000 else Nval=opt[['N_value']]
datagwas[,'N_gwas']<-Nval
headN<-'N_gwas'
}
datagwas$h2.gwas<-computedher(datagwas[,headbeta], datagwas[,headse], datagwas[,headaf],datagwas[,headN])
datagwas$z.gwas<-datagwas[,headbeta]/datagwas[,headse]
#datagwasf<-datagwas[datagwas[,headpval]<headpval,]

datalda1<-merge(datagwascat, datald, by.x=c(headchrcat,headbpcat), by.y=c("CHR_A", "BP_A"));names(datalda1)[names(datalda1)=="CHR_B"]<-headchr;names(datalda1)[names(datalda1)=="BP_B"]<-headbp;names(datalda1)[names(datalda1)=="SNP_B"]<-'rs_gwas';names(datalda1)[names(datalda1)=="SNP_A"]<-'rs_cat'
datalda2<-merge(datagwascat, datald, by.x=c(headchrcat,headbpcat), by.y=c("CHR_B", "BP_B"));names(datalda2)[names(datalda2)=="CHR_A"]<-headchr;names(datalda2)[names(datalda2)=="BP_A"]<-headbp;names(datalda2)[names(datalda2)=="SNP_A"]<-'rs_gwas';names(datalda2)[names(datalda2)=="SNP_B"]<-'rs_cat'

dataldallcat<-rbind(datalda1,datalda2)

if(Test){
rsclump<-"12:2523697";rscat="rs55935819";chrocat<-12;bpcat<-2521579
datagwas[datagwas[,headchr]==12 & datagwas[,headbp]==2523697,]
datagwascat[datagwascat[,headchrcat]==chrocat & datagwascat[,headbpcat]==bpcat,]
datagwascat[datagwascat[,'name']==rscat ,]

dataldallcat[dataldallcat$rs_gwas==rsclump,]
dataldallcat[dataldallcat$rs_cat==rsclump,]
dataldallcat[dataldallcat$rs_cat==rscat,]
datald[datald$SNP_B==rsclump | datald$SNP_A==rsclump,]
datagwas[datagwas$CHR=='18' & datagwas$BP==48133241,]
datald[datald$SNP_B==rscat | datald$SNP_A==rscat,]
#datagwascat[datagwascat$name=='rs745821',]
dataldallcat[dataldallcat$rs_cat==rscat,]
dataldallcat[dataldallcat$rs_gwas==rsclump,]
}


datalda1<-merge(dataldallcat, datagwas, by.x=c(headchr,headbp), by.y=c(headchr,headbp))
if(Test)datalda1[datalda1$name==rscat,]
#datald[datald$CHR_A==18 & datald$BP_A==48144571 ,]
#datalda1[datalda1$rs_cat=='rs745821',]
#datalda1[datalda1$rs_gwas=='18:48144571',]


write.table(datalda1, file=paste(opt[['out']],'_all.txt',sep=''), row.names=F, col.names=T,quote=F)




infocat=strsplit(opt[['info_gwascat']],split=';')[[1]]

datalda1$info_gwas<-paste(datalda1[,headchr],':',datalda1[,headbp],'-beta:',datalda1[,headbeta], ',se:',datalda1[,headse],',pval:',datalda1[,headpval],',R2:', datalda1[,'R2'],sep='')
datalda1$info_gwascat<-""
for(cat in infocat)datalda1$info_gwascat<-paste(datalda1$info_gwascat,cat,':',datalda1[,cat],',',sep='')
datagwassumm<-aggregate(as.formula(paste('info_gwas~',headbpcat, '+',headchrcat)), data=datalda1,function(x)paste(unique(x), collapse=';'))
datagwascatsumm<-aggregate(as.formula(paste('info_gwascat~',headbpcat, '+',headchrcat)), data=datalda1, function(x)paste(unique(x), collapse=';'))
datagwasminpval<-aggregate(as.formula(paste(headpval,'~',headbpcat, '+',headchrcat)), data=datalda1,min)
datagwasminpval<-merge(datagwasminpval,datalda1, by=c(headchrcat,headbpcat,headpval),all=F)


allresume<-merge(merge(datagwassumm,datagwascatsumm,all=T, by=c(headchrcat,headbpcat)),datagwasminpval, all=T, by=c(headchrcat,headbpcat))
names(allresume)[c(1,2)]<-c('chr_gwas', 'bp_gwas_cat')
write.csv(allresume, file=paste(opt[['out']],'_resume.csv',sep=''),row.names=F)

###
#head(datalda1)
#print(range(datalda1$R2))
#cat(opt[['min_pvalue']])
datalda1sig<-datalda1[datalda1[,headpval]<opt[['min_pvalue']],]
#datalda1sig<-datalda1
#if(nrow(datalda1sig)>0){
#datalda1sig$info_gwas<-paste(datalda1sig[,headchr],':',datalda1sig[,headbp],'-beta:',datalda1sig[,headbeta], ',se:',datalda1sig[,headse],',pval:',datalda1sig[,headpval],',R2:', datalda1sig[,'R2'],sep='')
#}else{
#}
#datalda1sig$info_gwascat<-""
#for(cat in infocat)datalda1sig$info_gwascat<-paste(datalda1sig$info_gwascat,cat,':',datalda1sig[,cat],',',sep='')
if(nrow(datalda1sig)>0){
datagwasminpval<-aggregate(as.formula(paste(headpval,'~',headbpcat, '+',headchrcat)), data=datalda1sig,min)
datagwassumm<-aggregate(as.formula(paste('info_gwas~',headbpcat, '+',headchrcat)), data=datalda1sig,function(x)paste(unique(x), collapse=';'))
datagwascatsumm<-aggregate(as.formula(paste('info_gwascat~',headbpcat, '+',headchrcat)), data=datalda1sig, function(x)paste(unique(x), collapse=';'))
}else{
datagwasminpval<-datalda1sig[F,c(headbpcat,headchrcat,headpval)]
datagwassumm<-datalda1sig[F,c('info_gwas',headchrcat,headbpcat)]
datagwascatsumm<-datalda1sig[F,c('info_gwascat',headchrcat,headbpcat)]
}


allresume<-merge(merge(datagwassumm,datagwascatsumm,all=T, by=c(headchrcat,headbpcat)), datagwasminpval,all=T, by=c(headchrcat,headbpcat))
allresume<-allresume[allresume[,headpval]<opt[['min_pvalue']],]
names(allresume)[c(1,2)]<-c('chr_gwas', 'bp_gwas_cat')
#minpval<-aggregate(as.formula(paste(headpval,'~',headbpcat, '+',headchrcat)), data=datalda1sig,min)
#names(minpval)[3]<-"min_pvalgwas"
#allresume<-merge(allresume,minpval,by=c(1,2))

## write sig
write.csv(allresume, file=paste(opt[['out']],'_resumesig.csv',sep=''),row.names=F)


