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

#--chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#opt=list(gwascat='Diastolic_AwigenLD_all.csv',gwas='Diastolic_AwigenLD_pos.init',chr_gwas='CHR',ps_gwas='BP',a1_gwas='ALLELE1',a2_gwas='ALLELE0',beta_gwas='BETA',se_gwas='SE',af_gwas='A1FREQ',chr_gwascat='chrom',bp_gwascat='chromEnd',p_gwas='P_BOLT_LMM',ps_gwascat='chromEnd',chr_gwascat='chrom',out='Diastolic_AwigenLD_ld',ld_file='Diastolic_AwigenLD_ld.ld',min_pvalue='0.001',min_r2=0.2,info_gwascat="pubMedID;author;trait;initSample")


headse=opt[['se_gwas']];headbp=opt[['ps_gwas']];headchr=opt[['chr_gwas']];headbeta=opt[['beta_gwas']];heada1=opt[['a1_gwas']];heada2=opt[['a2_gwas']];headpval=opt[['p_gwas']];headaf<-opt[['af_gwas']];headbeta=opt[['beta_gwas']]
headchrcat=opt[['chr_gwascat']];headbpcat=opt[['ps_gwascat']];heada1catrs<-"riskAllele";headzcat="z.cat";headafcat<-'risk.allele.af';heada1cat<-'risk.allele.cat'
outhead=opt[['out']]


datagwascat=read.csv(opt[['gwascat']])
datagwascat[,heada1cat]<-sapply(strsplit(as.character(datagwascat[,heada1catrs]),split='-'),function(x)x[2])
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
datagwasf<-datagwas[datagwas[,headpval]<headpval,]

datalda1<-merge(datagwascat, datald, by.x=c(headchrcat,headbpcat), by.y=c("CHR_A", "BP_A"));names(datalda1)[names(datalda1)=="CHR_B"]<-headchr;names(datalda1)[names(datalda1)=="BP_B"]<-headbp;names(datalda1)[names(datalda1)=="SNP_B"]<-'rs_gwas';names(datalda1)[names(datalda1)=="SNP_A"]<-'rs_cat'
datalda2<-merge(datagwascat, datald, by.x=c(headchrcat,headbpcat), by.y=c("CHR_B", "BP_B"));names(datalda2)[names(datalda2)=="CHR_A"]<-headchr;names(datalda2)[names(datalda2)=="BP_A"]<-headbp;names(datalda2)[names(datalda2)=="SNP_A"]<-'rs_gwas';names(datalda2)[names(datalda2)=="SNP_B"]<-'rs_cat'

dataldallcat<-rbind(datalda1,datalda2)
#rs61148001	NApubMedID:28135244,author:Warren HR,trait:Diastolic blood pressure,initSample:140,886 European ancestry individuals,	rs745821	18:48111983-beta:0.686852,se:0.245673,pval:0.0052;18:48122758-beta:0.783012,se:0.274825,pval:0.0044;18:48142854-beta:0.544822,se:0.255059,pval:0.033;18:48130202-beta:0.784569,se:0.275135,pval:0.0044;18:48109325-beta:-0.258781,se:0.220929,pval:0.24;18:48116645-beta:0.681003,se:0.2466,pval:0.0058;18:48146048-beta:0.4173,se:0.239155,pval:0.081;18:48108763-beta:-0.262851,se:0.222183,pval:0.24;18:48124420-beta:0.783787,se:0.274543,pval:0.0043;18:48123761-beta:-0.453955,se:0.21873,pval:0.038;18:48144571-beta:1.04508,se:0.282915,pval:0.00022;18:48142784-beta:0.525746,se:0.254897,pval:0.039;18:48124763-beta:0.799733,se:0.27479,pval:0.0036;18:48111045-beta:0.687581,se:0.245567,pval:0.0051;18:48133241-beta:1.05306,se:0.278833,pval:0.00016;18:48141710-beta:1.04796,se:0.282118,pval:2e-04;18:48111441-beta:0.783753,se:0.274781,pval:0.0043;18:48130532-beta:0.516826,se:0.2316,pval:0.026;18:48138375-beta:0.715874,se:0.265017,pval:0.0069;18:48140734-beta:0.767078,se:0.265891,pval:0.0039;18:48147127-beta:0.524988,se:0.255126,pval:0.04;18:48140238-beta:0.533419,se:0.255851,pval:0.037;18:48132646-beta:0.887984,se:0.270482,pval:0.001;18:48115431-beta:0.672788,se:0.245552,pval:0.0062;18:48111941-beta:0.779641,se:0.274712,pval:0.0045;18:48125899-beta:0.794057,se:0.274877,pval:0.0039	rs1370479;rs73959979;rs745821;rs73959982;rs4939637;rs73959977;rs3867257;rs2969972;rs73959980;rs930362;rs17742138;rs745823;rs17742008;rs10432161;rs61148001;rs58693787;rs10432162;rs10502909;rs58000211;rs11082833;rs1025686;rs4599004;rs1437649;rs11082832;rs1370480;rs8098724	18	48133241	25	3

datald[datald$SNP_B=='rs61148001' | datald$SNP_A=='rs61148001',]
datagwas[datagwas$CHR=='18' & datagwas$BP==48133241,]
datald[datald$SNP_B=='rs61148001' | datald$SNP_A=='rs61148001',]
#datagwascat[datagwascat$name=='rs745821',]
#dataldallcat[dataldallcat$rs_cat=='rs745821',]

#> datald[datald$CHR_A==18 & datald$BP_A==48144571 ,]
#   CHR_A     BP_A      SNP_A CHR_B     BP_B      SNP_B       R2
#1:    18 48144571 rs17742138    18 48146048  rs3867257 0.541578
#2:    18 48144571 rs17742138    18 48144571 rs17742138 1.000000
#> datald[datald$CHR_B==18 & datald$BP_B==48144571 ,]
#   CHR_A     BP_A      SNP_A CHR_B     BP_B      SNP_B R2
#1:    18 48144571 rs17742138    18 48144571 rs17742138  1


datalda1<-merge(dataldallcat, datagwas, by.x=c(headchrcat,headbpcat), by.y=c(headchr,headbp))
#datald[datald$CHR_A==18 & datald$BP_A==48144571 ,]
#datalda1[datalda1$rs_cat=='rs745821',]
#datalda1[datalda1$rs_gwas=='18:48144571',]


write.table(datalda1, file=paste(opt[['out']],'_all.txt',sep=''), row.names=F, col.names=T,quote=F)




infocat=strsplit(opt[['info_gwascat']],split=';')[[1]]

datalda1$info_gwas<-paste(datalda1[,headchr],':',datalda1[,headbp],'-beta:',datalda1[,headbeta], ',se:',datalda1[,headse],',pval:',datalda1[,headpval],sep='')
datalda1$info_gwascat<-""
for(cat in infocat)datalda1$info_gwascat<-paste(datalda1$info_gwascat,cat,':',datalda1[,cat],',',sep='')
datagwassumm<-aggregate(as.formula(paste('info_gwas~',headbpcat, '+',headchrcat)), data=datalda1,function(x)paste(unique(x), collapse=';'))
datagwascatsumm<-aggregate(as.formula(paste('info_gwascat~',headbpcat, '+',headchrcat)), data=datalda1, function(x)paste(unique(x), collapse=';'))

allresume<-merge(datagwassumm,datagwascatsumm,all=T, by=c(headchrcat,headbpcat))
names(allresume)[c(1,2)]<-c('chr_gwas', 'bp_gwas_cat')
write.csv(allresume, file=paste(opt[['out']],'_resume.csv',sep=''),row.names=F)

###
head(datalda1)
print(range(datalda1$R2))
cat(opt[['min_pvalue']])
datalda1sig<-datalda1[datalda1[,headpval]<opt[['min_pvalue']],]
#datalda1sig<-datalda1
datalda1sig$info_gwas<-paste(datalda1sig[,headchr],':',datalda1sig[,headbp],'-beta:',datalda1sig[,headbeta], ',se:',datalda1sig[,headse],',pval:',datalda1sig[,headpval],',R2:', datalda1sig[,'R2'],sep='')
datalda1sig$info_gwascat<-""
for(cat in infocat)datalda1sig$info_gwascat<-paste(datalda1sig$info_gwascat,cat,':',datalda1sig[,cat],',',sep='')
datagwasminpval<-aggregate(as.formula(paste(headpval,'~',headbpcat, '+',headchrcat)), data=datalda1sig,min)

datagwassumm<-aggregate(as.formula(paste('info_gwas~',headbpcat, '+',headchrcat)), data=datalda1sig,function(x)paste(unique(x), collapse=';'))
datagwascatsumm<-aggregate(as.formula(paste('info_gwascat~',headbpcat, '+',headchrcat)), data=datalda1sig, function(x)paste(unique(x), collapse=';'))

allresume<-merge(merge(datagwassumm,datagwascatsumm,all=T, by=c(headchrcat,headbpcat)), datagwasminpval,all=T, by=c(headchrcat,headbpcat))
allresume<-allresume[allresume[,headpval]<opt[['min_pvalue']],]
names(allresume)[c(1,2)]<-c('chr_gwas', 'bp_gwas_cat')
minpval<-aggregate(as.formula(paste(headpval,'~',headbpcat, '+',headchrcat)), data=datalda1sig,min)
#names(minpval)[3]<-"min_pvalgwas"
#allresume<-merge(allresume,minpval,by=c(1,2))

## write sig
write.csv(allresume, file=paste(opt[['out']],'_resumesig.csv',sep=''),row.names=F)


