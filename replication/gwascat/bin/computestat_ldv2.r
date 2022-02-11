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
  make_option(c( "--ldblock_file"), type="character", default=NULL,
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
opt = parse_args(opt_parser)
if(Test)opt = parse_args(opt_parser, strsplit('--gwascat DBP_AwigenLD_all.csv --gwas DBP_AwigenLD_range.init --chr_gwas CHR --ps_gwas BP --a1_gwas ALLELE1 --a2_gwas ALLELE0 --beta_gwas BETA --se_gwas SE --chr_gwascat chrom --bp_gwascat chromEnd --p_gwas P_BOLT_LMM --ps_gwascat chromEnd --chr_gwascat chrom --out DBP_AwigenLD_wind --ldblock_file tmpout --min_pvalue 0.01 --min_r2 0.5 --info_gwascat "pubMedID;author;trait;initSample', split=' ')[[1]])


#/home/jeantristan/Travail/git/h3agwas/replication/gwascat/bin/computestat_ldv2.r  --gwascat meanMaxcIMT_eurld_all.csv --gwas meanMaxcIMT_eurld_range.init --chr_gwas CHR --ps_gwas BP --a1_gwas ALLELE1 --a2_gwas ALLELE0  --beta_gwas BETA --se_gwas SE  --chr_gwascat chrom --bp_gwascat chromEnd --p_gwas P_BOLT_LMM --ps_gwascat chromEnd --chr_gwascat chrom --out meanMaxcIMT_eurld_noext --ldblock_file meanMaxcIMT_eurld_ldinfo_ld.out --min_pvalue 1 --min_r2  0.5 --info_gwascat "pubMedID;author;trait;initSample"
#if(Test)opt=list(gwascat='meanMaxcIMT_eurld_all.csv',gwas='meanMaxcIMT_eurld_range.init',chr_gwas='CHR',ps_gwas='BP',a1_gwas='ALLELE1',a2_gwas='ALLELE0',beta_gwas='BETA',se_gwas='SE',af_gwas='A1FREQ',chr_gwascat='chrom',bp_gwascat='chromEnd',p_gwas='P_BOLT_LMM',ps_gwascat='chromEnd',chr_gwascat='chrom',out='out',ldblock_file='meanMaxcIMT_eurld_ldinfo_ld.out',min_pvalue='0.0001',min_r2=0.5,info_gwascat="pubMedID;author;trait;initSample")
#n/computestat_ldv2.r  --gwascat meanMaxcIMT_eurld_AwigenLD_all.csv --gwas meanMaxcIMT_eurld_AwigenLD_range.init --chr_gwas CHR --ps_gwas BP --a1_gwas ALLELE1 --a2_gwas ALLELE0  --beta_gwas BETA --se_gwas SE  --chr_gwascat chrom --bp_gwascat chromEnd --p_gwas P_BOLT_LMM --ps_gwascat chromEnd --chr_gwascat chrom --out meanMaxcIMT_eurld_AwigenLD_ld2 --ldblock_file meanMaxcIMT_eurld_AwigenLD_ld2.tmp_pos --min_pvalue 1.0E-4 --min_r2  0.5 --info_gwascat "pubMedID;author;trait;initSample"



headse=opt[['se_gwas']];headbp=opt[['ps_gwas']];headchr=opt[['chr_gwas']];headbeta=opt[['beta_gwas']];heada1=opt[['a1_gwas']];heada2=opt[['a2_gwas']];headpval=opt[['p_gwas']];headaf<-opt[['af_gwas']];headbeta=opt[['beta_gwas']]
headchrcat=opt[['chr_gwascat']];headbpcat=opt[['ps_gwascat']];heada1catrs<-"riskAllele";headzcat="z.cat";headafcat<-'risk.allele.af';heada1cat<-'risk.allele.cat'
outhead=opt[['out']]

version<-2

datagwascat=read.csv(opt[['gwascat']])
#datagwascat[,heada1cat]<-sapply(strsplit(as.character(datagwascat[,heada1catrs]),split='-'),function(x)x[2])
datagwas<-fread(opt[['gwas']], header=T)
if(!is.null(headaf))checkhead(headaf, datagwas,'af')
checkhead(headpval, datagwas,'pval');checkhead(headse, datagwas,'se');checkhead(headbp, datagwas,'bp');checkhead(headchr, datagwas, 'chr');checkhead(headbeta, datagwas, 'beta')

checkhead(headbpcat,datagwascat,'bp cat');checkhead(headchrcat,datagwascat,'chro cat');

# CHR_A         BP_A         SNP_A  CHR_B         BP_B         SNP_B           R2 
data_resumld<-fread(opt[['ldblock_file']])
if(version==1)names(data_resumld)<-c('block', 'CHR', 'BP', 'TYPE') else names(data_resumld)<-c('block', 'CHR', 'BP')
data_resumld$block<-paste(data_resumld$CHR,data_resumld$block, sep='_')

#tmpa<-unique(datald[,c('CHR_A','BP_A','SNP_A')]);names(tmpa)<-c('CHR','BP','SNP');
#tmpb<-unique(datald[,c('CHR_B','BP_B','SNP_B')]);names(tmpb)<-c('CHR','BP','SNP');
#tmp<-rbind(tmpa,tmpb)
#tmpall<-unique(rbind(tmpa,tmpb))
#tmpall<-tmpall[,c('CHR','BP','SNP','CHR','BP','SNP')]
#tmpall$R2<-1
#names(tmpall)<-names(datald)
#datald<-rbind(datald,tmpall)

#   CHR_A     BP_A     SNP_A CHR_B     BP_B      SNP_B       R2
#1:    18 48132646 rs1437649    18 48133241 rs61148001 0.896039





headN<-opt[['N_value']]
if(is.null(opt[['N_gwas']])){
if(is.null(opt[['N_value']]))Nval<-10000 else Nval=opt[['N_value']]
datagwas[['N_gwas']]<-Nval
headN<-'N_gwas'
}
if(!is.null(headaf)){
 datagwas$h2.gwas<-computedher(datagwas[[headbeta]], datagwas[[headse]], datagwas[[headaf]],datagwas[[headN]])
 datagwas$z.gwas<-datagwas[[headbeta]]/datagwas[[headse]]
}else{
datagwas$h2.gwas<-NA
datagwas$z.gwas<-NA

}

#datalda1<-merge(datagwascat, datald, by.x=c(headchrcat,headbpcat), by.y=c("CHR_A", "BP_A"));names(datalda1)[names(datalda1)=="CHR_B"]<-headchr;names(datalda1)[names(datalda1)=="BP_B"]<-headbp;names(datalda1)[names(datalda1)=="SNP_B"]<-'rs_gwas';names(datalda1)[names(datalda1)=="SNP_A"]<-'rs_cat'
#datalda2<-merge(datagwascat, datald, by.x=c(headchrcat,headbpcat), by.y=c("CHR_B", "BP_B"));names(datalda2)[names(datalda2)=="CHR_A"]<-headchr;names(datalda2)[names(datalda2)=="BP_A"]<-headbp;names(datalda2)[names(datalda2)=="SNP_A"]<-'rs_gwas';names(datalda2)[names(datalda2)=="SNP_B"]<-'rs_cat'

#dataldallcat<-rbind(datalda1,datalda2)

dataresall<-merge(data_resumld,datagwascat, all.x=T, by.x=c('CHR', 'BP'), by.y=c(headchrcat,headbpcat),allow.cartesian=TRUE)
dataresall<-as.data.frame(merge(dataresall, datagwas, all.x=T, by.x=c('CHR', 'BP'), by.y=c(headchr, headbp), allow.cartesian=TRUE))


write.table(dataresall, file=paste(opt[['out']],'_all.txt',sep=''), row.names=F, col.names=T,quote=F, sep='\t')

infocatI=strsplit(opt[['info_gwascat']],split=';')[[1]]

dataresall$info_gwas<-paste(dataresall[,'CHR'],':',dataresall[,'BP'],'-beta:',dataresall[,headbeta], ',se:',dataresall[,headse],',pval:',dataresall[,headpval])
dataresall$info_gwascat<-""
balise<-!apply(dataresall[,infocatI], 1, function(x)all(is.na(x)))
infocat<-c('CHR','BP',infocatI)
for(cat in infocat)dataresall$info_gwascat[balise]<-paste(dataresall$info_gwascat[balise],cat,':',dataresall[balise,cat],',',sep='')
write.csv(dataresall, file=paste(opt[['out']],'_all.csv',sep=''),row.names=F)

infocatdata<-aggregate(info_gwascat~block, data=dataresall[,],function(x)paste(unique(x), collapse=';'))
infodata<-aggregate(info_gwas~block, data=dataresall,function(x)paste(unique(x), collapse=';'))
chro<-aggregate(as.formula(paste("CHR~block")), data=dataresall, unique)
bpmin<-aggregate(as.formula(paste("BP~block")), data=dataresall, min)
bpmax<-aggregate(as.formula(paste("BP~block")), data=dataresall, max)
#wh.minpvaldata<-aggregate(as.formula(paste('cbind(',headbp,',',headpval,')',"~block")), data=dataresall,function(x)print(x))

infobloc<-merge(merge(chro, bpmin, by='block',allow.cartesian=TRUE),bpmax, by='block',allow.cartesian=TRUE)

minpvaldata<-aggregate(as.formula(paste(headpval,"~block")), data=dataresall, min);names(minpvaldata)<-c('block', 'min_pval')
dataresalltmp<-merge(dataresall,minpvaldata, by='block')
minpvalbp<-dataresalltmp[dataresalltmp$min_pval==dataresalltmp[,headpval],]
minpos<-aggregate(BP~block ,data=minpvalbp,function(x)paste(x,collapse=';'))
names(minpos)<-c('block', 'best_pos')


print(head(infobloc))
names(infobloc)<-c('block', 'chro','min_bp', 'max_bp')


ndata<-aggregate(as.formula(paste(headpval,"~block")), data=dataresall, length)
names(ndata)<-c('block', 'n_total')

allmerge<-merge(merge(infobloc,merge(merge(merge(infocatdata,infodata,by='block',all=T),minpvaldata,by='block',all=T), ndata, by='block',all=T),by='block',all=T), minpos, by='block',all=T)
print(is.data.frame(allmerge))
write.csv(allmerge, file=paste(opt[['out']],'_allresume.csv',sep=''),row.names=F)

minpval<-as.numeric(opt[['min_pvalue']])
dataresallsig<-dataresall[!is.na(dataresall[,headpval]) & dataresall[,headpval]<minpval,]
if(nrow(dataresallsig)>0){
dataresallsig$info_gwas<-paste(dataresallsig[,'CHR'],':',dataresallsig[,'BP'],'-beta:',dataresallsig[,headbeta], ',se:',dataresallsig[,headse],',pval:',dataresallsig[,headpval])

chro<-aggregate(as.formula(paste("CHR~block")), data=dataresallsig, unique)

ndata<-aggregate(as.formula(paste(headpval,"~block")), data=dataresallsig, length)
names(ndata)<-c('block', 'n_total')

infodata<-aggregate(info_gwas~block, data=dataresallsig,function(x)paste(unique(x), collapse=';'))
minpvaldata<-aggregate(as.formula(paste(headpval,"~block")), data=dataresallsig, min)
ndataSig<-aggregate(as.formula(paste(headpval,"~block")), data=dataresallsig, length)

infobloc<-merge(merge(chro, bpmin, by='block'),bpmax, by='block')
names(infobloc)<-c('block', 'chro','min_bp', 'max_bp')

minpvaldata<-aggregate(as.formula(paste(headpval,"~block")), data=dataresallsig, min);names(minpvaldata)<-c('block', 'min_pval')
dataresalltmp<-merge(dataresallsig,minpvaldata, by='block')
minpvalbp<-dataresalltmp[dataresalltmp$min_pval==dataresalltmp[,headpval],]
minpos<-aggregate(BP~block ,data=minpvalbp,function(x)paste(x,collapse=';'))
names(minpos)<-c('block', 'best_pos')


names(ndataSig)<-c('block', 'n_sig')

allmergesig<-merge(merge(infobloc,merge(merge(merge(merge(infocatdata,infodata,by='block',all=F),minpvaldata,by='block',all=F), ndata, by='block',all=F), ndataSig, by='block', all=F),by='block',all=F), minpos,by='block',all=F)
write.csv(allmergesig, file=paste(opt[['out']],'_resumesig.csv',sep=''),row.names=F)
}else{
head=paste(c("block","chro","min_bp","max_bp","info_gwascat","info_gwas","min_pval","n_total","best_pos"), collapse=',')
writeLines(head, con=paste(opt[['out']],'_resumesig.csv',sep=''))

}




