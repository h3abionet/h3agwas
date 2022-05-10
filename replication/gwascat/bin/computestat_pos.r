#!/usr/bin/env Rscript
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

gopt<-function(x){
if(is.null(opt[[x]]))return(NULL)
gsub('-','.',opt[[x]])
}


plotfreq<-function(dataall,freq1, freq2,cex_pt,alpha_pt,xlab='GWAS frequencies', ylab='GWAS cat frequencies'){
nbcat=20
height=1
perctrans<-90
ht<-hist(dataall[,freq1], nbcat, plot=F)
ht$counts<-ht$counts/sum(ht$counts)*height
ht2<-hist(dataall[,freq2], nbcat, plot=F)
ht2$counts<-ht2$counts/sum(ht2$counts)*height
layout.matrix <- matrix(c(2, 1, 0, 3), nrow = 2, ncol = 2)
dataall<-na.omit(dataall[,c(freq1, freq2)])

layout(mat = layout.matrix,
       heights = c(max(ht$counts)*1.2, 1), # Heights of the two rows
       widths = c(2, max(ht2$counts)*2.4)) # Widths of the two columns
#layout.show(3)

#c(bottom, left, top, right)
par(mar=c(5, 4, 1, 1))
plot(c(0,1), c(0,1), type='n', cex=0.5, xlab=xlab, ylab=ylab, xlim=c(0,1), ylim=c(0,1), bty='n')
points(dataall[,freq1], dataall[,freq2], pch=22, ,col=t_col("blue",alpha_pt), bg=t_col("blue",alpha_pt),cex=cex_pt)
abline(a=0,b=1, lty=2 ,col=t_col('red'), lwd=2)
par(mar=c(0,4,1,0))
plot(ht, col=t_col("orange"), xlab="", ylab="", xlim=c(0,1),  xaxt='n',main="")
par(mar=c(5,0,1,1))
plot(c(0, max(ht2$counts)), c(0,1), type='n', cex=0.5, xlab='', ylab='', xlim=c(0,max(ht2$counts)), ylim=c(0,1), yaxt='n',bty='n')
rect(rep(0,length(ht2$counts)),(1:length(ht2$counts))/nbcat,
           ht2$counts,(1:length(ht2$counts))/nbcat-1/nbcat, col=t_col('red', perctrans), bg=t_col('red', perctrans))
}

plotZ<-function(dataall,Z1, Z2, cex_pt=0.5,alpha_pt=10,xlab='Z (GWAS cat)', ylab='Z (AWIGEN)'){
balise<-!is.na(dataall[,Z1]) & !is.na(dataall[,Z2])
dataall<-dataall[balise,]
r2<-cor(dataall[,Z1],dataall[,Z2], method='spearman')
r2abs<-cor(abs(dataall[,Z1]), abs(dataall[,Z2]), method='spearman')
plot(dataall[,Z1], dataall[,Z2], pch=22, cex=cex_pt,bg=t_col("blue", alpha_pt) ,col=t_col("blue",alpha_pt), xlab=xlab, ylab=ylab)
text(min(dataall[,Z1])-min(dataall[,Z1])*0.1,max(dataall[,Z2]), paste(paste("r2 :",round(r2,2)), paste("\n         r2 (abs) :",round(r2abs,2)), sep=""))
abline(h=0, col='red', lty=2)
abline(v=0, col='red', lty=2)
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
  make_option(c("--z_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--af_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--N_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--ps_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--chr_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--a1_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--p_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

checkhead<-function(head,Data, type){
if(length(which(head %in% names(Data)))==0){
print(names(Data))
print(paste('not found ', head,'for',type ,'data'))
q('no',2)
}
}

#--chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
Test=F
if(Test)opt=list(gwascat='meanMaxcIMT_eurld_all.csv',gwas='meanMaxcIMT_eurld_pos.init',chr_gwas='CHR',ps_gwas='BP',a1_gwas='ALLELE1',a2_gwas='ALLELE0' ,beta_gwas='BETA',se_gwas='SE',af_gwas='A1FREQ',chr_gwascat='chrom',bp_gwascat='chromEnd',p_gwas='P_BOLT_LMM',ps_gwascat='chromEnd',chr_gwascat='chrom',out='meanMaxcIMT_eurld_pos')


headse=gopt('se_gwas');headps=gopt('ps_gwas');headchr=gopt('chr_gwas');heada1=gopt('a1_gwas');heada2=gopt('a2_gwas');headpval=gopt('p_gwas');headaf<-gopt('af_gwas');headbeta=gopt('beta_gwas');headz=gopt('z_gwas')
headchrcat=gopt('chr_gwascat');headbpcat=gopt('ps_gwascat');heada1catrs<-gopt('a1_gwascat');headzcat="z.cat";headafcat<-'risk.allele.af';heada1cat<-'risk.allele.cat'

outhead=opt[['out']]


datagwascat=read.csv(opt[['gwascat']])
#datagwascat[,heada1cat]<-sapply(strsplit(as.character(datagwascat[,heada1catrs]),split='-'),function(x)x[2])
datagwas<-read.table(opt[['gwas']], header=T)
checkhead(headpval, datagwas,'pval');checkhead(headps, datagwas,'bp');checkhead(headchr, datagwas, 'chr');

if(!is.null(headaf)){
checkhead(headaf, datagwas,'af')
}
if(!is.null(headbeta) & !is.null(headse)){
print('beta, se not null')
checkhead(headbeta, datagwas, 'beta');checkhead(headse, datagwas,'se')
}
if(!is.null(headz)){
checkhead(headz, datagwas, 'z_gwas')
datagwas$z.gwas<-datagwas[[headz]]
}

checkhead(headbpcat,datagwascat,'bp cat');checkhead(headchrcat,datagwascat,'chro cat');


print(head(datagwas,1))
datagwas$p.adjust.fdr<-p.adjust(datagwas[,headpval],'fdr')
datagwas$p.adjust.bonf<-p.adjust(datagwas[,headpval],'bonferroni')
headN<-opt[['N_value']]
if(is.null(opt[['N_gwas']])){
if(is.null(opt[['N_value']]))Nval<-10000 else Nval=opt[['N_value']]
datagwas[,'N_gwas']<-Nval
headN<-'N_gwas'
}

if(!is.null(headbeta))datagwas$z.gwas<-datagwas[,headbeta]/datagwas[,headse]
if(!is.null(headaf) & !is.null(headbeta)){
datagwas$h2.gwas<-computedher(datagwas[,headbeta], datagwas[,headse], datagwas[,headaf],datagwas[,headN])
}else{
datagwas$h2.gwas<-NA
}


MergeAll<-merge(datagwas, datagwascat,by.x=c(headchr,headps), by.y=c(headchrcat,headbpcat))

MergeAll[,heada1cat]<-as.character(MergeAll[,heada1cat])
MergeAll[,heada1]<-as.character(MergeAll[,heada1])

QC<-!is.na(MergeAll[,heada1]) & !is.na(MergeAll[,heada2]) & !is.na(MergeAll[,heada1cat]) & !is.na(MergeAll[,headafcat]) & (MergeAll[,heada1] == MergeAll[,heada1cat] | MergeAll[,heada2] == MergeAll[,heada1cat]) 
BaliseChange<-QC & MergeAll[,heada1]!=MergeAll[,heada1cat]

MergeAll[,'af_gwas_a1cat']<-NA
MergeAll[BaliseChange,'af_gwas_a1cat']<- 1 - MergeAll[BaliseChange,headafcat]
MergeAll[!BaliseChange,'af_gwas_a1cat']<-  MergeAll[!BaliseChange,headafcat]


pdf(paste(outhead,'_cmpfrequencies.pdf', sep=''))
if(!is.null(headaf)){
QC<-!is.na(MergeAll[,headaf]) & !is.na(MergeAll[,'af_gwas_a1cat'])
print(MergeAll[,headaf])
print(MergeAll[,'af_gwas_a1cat'])
plotfreq(MergeAll[QC,],headaf, 'af_gwas_a1cat',cex_pt=1.5,alpha_pt=0.15,xlab='GWAS', ylab='GWAS Catalog')
}else{
plot(1:10, 1:10, type='n')
text(2,5, 'no frequencie to plot')
}
dev.off()
QC<-!is.na(MergeAll[,heada1]) & !is.na(MergeAll[,heada2]) & !is.na(MergeAll[,heada1cat]) & !is.na(MergeAll[,headzcat]) & (MergeAll[,heada1] == MergeAll[,heada1cat] | MergeAll[,heada2] == MergeAll[,heada1cat])

BaliseChange<-QC & MergeAll[,heada1]!=MergeAll[,heada1cat]
MergeAll[,'z_gwas_a1cat']<-NA
MergeAll[BaliseChange,'z_gwas_a1cat']<- - MergeAll[BaliseChange,headzcat]
MergeAll[!BaliseChange,'z_gwas_a1cat']<-MergeAll[!BaliseChange,headzcat]


pdf(paste(outhead,'_cmpz.pdf', sep=''))
if(!is.null(headaf)){
print(head(MergeAll))
QC<-!is.na(MergeAll[,'z.gwas']) & !is.na(MergeAll[,'z_gwas_a1cat'])
plotZ(MergeAll[QC,],'z.gwas','z_gwas_a1cat', cex_pt=1.5,alpha_pt=0.15,ylab='GWAS catalog', xlab='GWAS')
}else{
plot(1:10, 1:10, type='n')
text(2,5, 'no frequencie to compute z ')
}
dev.off()

pdf(paste(outhead,'_qq.pdf', sep=''))
qq(MergeAll[,headpval])
dev.off()



write.csv(MergeAll, row.names=F,file=paste(opt[['out']], '.csv',sep=''))
