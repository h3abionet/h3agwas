#!/usr/bin/env Rscript
library("optparse")
library("qqman")
gopt<-function(x){
gsub('-','.',opt[[x]])
}
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


plotfreq<-function(dataall,freq1, freq2,cex_pt,alpha_pt,xlab='GWAS frequencies', ylab='GWAS cat frequencies'){
nbcat=20
height=1
perctrans<-90
ht<-hist(dataall[,freq1], nbcat, plot=F)
ht$counts<-ht$counts/sum(ht$counts)*height
ht2<-hist(dataall[,freq2], nbcat, plot=F)
ht2$counts<-ht2$counts/sum(ht2$counts)*height
layout.matrix <- matrix(c(2, 1, 0, 3), nrow = 2, ncol = 2)

layout(mat = layout.matrix,
       heights = c(max(ht$counts)*1.2, 1), # Heights of the two rows
       widths = c(2, max(ht2$counts)*2.4)) # Widths of the two columns
#layout.show(3)

#c(bottom, left, top, right)
par(mar=c(5, 4, 1, 1))
plot(c(0,1), c(0,1), type='n', cex=0.5, xlab=xlab, ylab=ylab, xlim=c(0,1), ylim=c(0,1), bty='n')
points(dataall[,freq1], dataall[,freq2], pch=22, ,col=t_col("blue",95), bg=t_col("blue",95),cex=cex_pt)
abline(a=0,b=1, lty=2 ,col=t_col('red'), lwd=2)
par(mar=c(0,4,1,0))
plot(ht, col=t_col("orange"), xlab="", ylab="", xlim=c(0,1),  xaxt='n',main="")
par(mar=c(5,0,1,1))
plot(c(0, max(ht2$counts)), c(0,1), type='n', cex=0.5, xlab='', ylab='', xlim=c(0,max(ht2$counts)), ylim=c(0,1), yaxt='n',bty='n')
rect(rep(0,length(ht2$counts)),(1:length(ht2$counts))/nbcat,
           ht2$counts,(1:length(ht2$counts))/nbcat-1/nbcat, col=t_col('red', perctrans), bg=t_col('red', perctrans))
}



plotZ<-function(dataall,Z1, Z2, xlab='Z (GWAS cat)', ylab='Z (AWIGEN)'){
datagwas$z.gwas<-NA
r2<-cor(dataall[,Z1],dataall[,Z2], method='spearman')
r2abs<-cor(abs(dataall[,Z1]), abs(dataall[,Z2]), method='spearman')
plot(dataall[,Z1], dataall[,Z2], pch=22, cex=0.5,bg=t_col("blue") ,col=t_col("blue"), xlab=xlab, ylab=ylab)
text(min(dataall[,Z1])-min(dataall[,Z1])*0.1,max(dataall[,Z2]), paste(paste("r2 :",round(r2,2)), paste("\n         r2 (abs) :",round(r2abs,2)), sep=""))
abline(h=0, col='red', lty=2)
abline(v=0, col='red', lty=2)
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
  make_option(c("--af_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--N_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--info_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--ps_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--chr_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--p_gwas"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--wind"), type="numeric", default=50,
              help="dataset file name", metavar="character"),
  make_option(c("--min_pval"), type="numeric", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--a1_gwascat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--merge_wind"), type="character", default=1,
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
#opt=list(gwascat='out_all.csv',gwas='out_range.init',chr_gwas='chr',ps_gwas='ps',a1_gwas='allele1',a2_gwas='allele0',beta_gwas='beta',se_gwas='se',af_gwas='af',chr_gwascat='chrom',bp_gwascat='ps',p_gwas='p_wald',ps_gwascat='chromEnd',chr_gwascat='chrom',out='out_pos', wind=25, min_pval=0.01,info_gwascat="pubMedID;author;trait;initSample")


headse=gopt('se_gwas');headbp=gopt('ps_gwas');headchr=gopt('chr_gwas');heada1=gopt('a1_gwas');heada2=gopt('a2_gwas');headpval=gopt('p_gwas');headaf<-gopt('af_gwas')
headchrcat=gopt('chr_gwascat');headbpcat=gopt('ps_gwascat');heada1catrs<-gopt('a1_gwascat');headzcat="z.cat";headafcat<-'risk.allele.af';heada1cat<-'risk.allele.cat'
outhead=opt[['out']]
wind=opt[['wind']]*1000


datagwascat=read.csv(opt[['gwascat']])
checkhead(headbpcat,datagwascat,'bp cat');checkhead(headchrcat,datagwascat,'chro cat')
#checkhead(opt[['wind']],datagwascat,'chro cat');

#datagwascat[,heada1cat]<-sapply(strsplit(as.character(datagwascat[,heada1catrs]),split='-'),function(x)x[2])
datagwascat$begin<-datagwascat[,headbpcat]-wind
datagwascat$end<-datagwascat[,headbpcat]+wind
datagwas<-read.table(opt[['gwas']], header=T)
checkhead(headpval, datagwas,'pval');checkhead(headbp, datagwas,'bp');checkhead(headchr, datagwas, 'chr')



#datagwas$p.adjust.fdr<-p.adjust(datagwas[,headpval],'fdr')
#datagwas$p.adjust.bonf<-p.adjust(datagwas[,headpval],'bonferroni')
headN<-opt[['N_value']]
if(is.null(opt[['N_gwas']])){
if(is.null(opt[['N_value']]))Nval<-10000 else Nval=opt[['N_value']]
datagwas[,'N_gwas']<-Nval
headN<-'N_gwas'
}
baliseh2=F

##  merge windows
if(opt[['merge_wind']]==1){
NumWind<-1
Cmt<-1
for(chro in unique(datagwascat[,headchrcat])){
datagwascatchr<-datagwascat[datagwascat[,headchrcat]==chro,]
datagwascatchr$wind_num<- -1
datagwascatchr<-datagwascatchr[order(datagwascatchr$begin),]
datagwascatchr$wind_num[1]<-NumWind
if(nrow(datagwascatchr)>1)
for(wind in 2:nrow(datagwascatchr)){
if(datagwascatchr$begin[wind]<datagwascatchr$end[wind-1]){
NumWind<-NumWind
}else{
NumWind<-NumWind+1
}
datagwascatchr$wind_num[wind]<-NumWind
}
if(Cmt==1)datagwascatf<-datagwascatchr
else datagwascatf<-rbind(datagwascatf ,datagwascatchr)
Cmt<-Cmt+1
}
}else{
datagwascatf<-datagwascat
datagwascatf$wind_num<-1:nrow(datagwascatf)

}
MinWind<-aggregate(as.formula(paste('begin~wind_num+',headchrcat,sep='')), data=datagwascatf, FUN=min)
MaxWind<-aggregate(as.formula(paste('end~wind_num+',headchrcat,sep='')), data=datagwascatf, FUN=max)

AllWind<-merge(MinWind,MaxWind)


AllWind[,headchrcat]<-as.character(AllWind[,headchrcat])
datagwas[,headchr]<-as.character(datagwas[,headchr])
datagwas$wind_num<-NA
for(Cmtwindcat in 1:nrow(AllWind)){
windcat=AllWind[Cmtwindcat,]
balise=datagwas[,headchr]==windcat[,headchrcat] & datagwas[,headbp] >=windcat$begin & datagwas[,headbp] <=windcat$end 
datagwas$wind_num[balise]<-windcat$wind_num
windcat$nbpos<-length(which(balise))
windcat$nbpos_sig<-length(which(balise & datagwas[,headpval]<opt[['min_pval']]))
PosMostSig<-which(min(datagwas[balise,headpval])==datagwas[,headpval] & balise)[1]
windcat$min_bp_gwas<-datagwas[PosMostSig,headbp]
if(Cmtwindcat==1)newwindcat<-windcat
else newwindcat <-rbind(newwindcat,windcat)
}

Allwingwasca<-merge(datagwascatf[,!(names(datagwascatf) %in% c("begin","end"))],merge(newwindcat,datagwas[!is.na(datagwas$wind_num) ,],by='wind_num'),by=c('wind_num'))

#write.csv(Allwingwasca,file=paste(opt[['out']],'_resume_allwind.csv',sep=''))
best_windcat<-merge(datagwascatf[,!(names(datagwascatf) %in% c("begin","end"))],merge(newwindcat,datagwas[!is.na(datagwas$wind_num) ,],by.x=c(headchrcat,'min_bp_gwas','wind_num'), by.y=c(headchr, headbp,'wind_num')),by=c('wind_num'))
write.csv(best_windcat,file=paste(opt[['out']],'_bestgwaswind.csv',sep=''))

infocat=strsplit(opt[['info_gwascat']],split=';')[[1]]
balise<-TRUE
best_windcat$info_gwascat=""
for(cat in infocat)best_windcat$info_gwascat[balise]<-paste(best_windcat$info_gwascat[balise],cat,':',best_windcat[balise,cat],',',sep='')

tmpinfocat<-aggregate(info_gwascat~wind_num, data=best_windcat,FUN=paste, collapse=';')
best_windcat_2<-merge(newwindcat,datagwas[!is.na(datagwas$wind_num) ,],by.x=c(headchrcat,'min_bp_gwas','wind_num'), by.y=c(headchr, headbp,'wind_num'))
resumebest<-merge(best_windcat_2,tmpinfocat, by='wind_num')

write.csv(resumebest,file=paste(opt[['out']],'_resume_bestgwaswind.csv',sep=''))




#Allwingwasca<-merge(datagwascatf[,names(datagwascatf)!=c("begin","end")],merge(AllWind,newwindcat, by='wind_num'),by=c('wind_num'))


