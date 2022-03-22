#!/usr/bin/env Rscript
library("optparse")


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
  make_option(c("--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--pheno"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--file_pheno"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--wind"), type="numeric", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--chro"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--chro_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--bp_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--pheno_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--beta_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--ci_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--p_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--n_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--freq_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--rs_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--riskall_head"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--format"), type="character", default="USCS", 
              help="dataset file name", metavar="character"),
  make_option(c("--typeformat"), type="character", default="csv",
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

balisepheno=F
if(!is.null(opt[['pheno']])){
lisp<-opt[['pheno']]
balisepheno=T
}else{
if(!is.null(opt[['file_pheno']])){
lisp<-readLines(opt[['file_pheno']])
balisepheno=T
}
}

lisc<-unlist(strsplit(opt[['chro']],','))
filegc<-opt[['file']]


#lisp="Type 2 diabetes";lisc=c("22","21")
format=opt[['format']]
if(format=="USCS"){
Data<-read.csv(opt[['file']], header=F, sep='\t')
chrohead="chrom";phenhead="trait";poshead="chromEnd";OrBetaHead="orOrBeta";headCi95="ci95";pValueHead="pValue";nvalueHead="initSample";freqHead="riskAlFreq";rsHead="name";riskall="riskAllele"
names(Data)<-c("bin", "chrom", "chromStart", "chromEnd", "name", "pubMedID", "author", "pubDate", "journal", "title", "trait", "initSample","replSample","region", "genes", "riskAllele", "riskAlFreq", "pValue", "pValueDesc", "orOrBeta", "ci95","platform", "cnv")
}else{
if(length(grep('.csv$', filegc))>0)Data<-read.csv(opt[['file']], header=T) else  Data<-read.csv(opt[['file']], header=T, sep="\t")
chrohead=opt[['chro_head']];phenhead=opt[['pheno_head']];poshead=opt[['bp_head']];OrBetaHead=opt[['beta_head']];headCi95=opt[['ci_head']];pValueHead=opt[['p_head']];nvalueHead=opt[['n_head']];freqHead=opt[['freq_head']];rsHead=opt[['rs_head']];riskall=opt[['riskall_head']]
}
Data[,chrohead]<-gsub('chr', '',as.character(Data[,chrohead]))
if(!is.null(opt[['chro']]))baliselistc<-Data[,chrohead] %in% lisc else baliselistc=T
if(balisepheno)balisepheno<-tolower(as.character(Data[,phenhead])) %in% tolower(lisp) else balisepheno=T
Data2<-Data[baliselistc & balisepheno,]
Data2[,'risk.allele.cat']<-sapply(strsplit(as.character(Data2[,riskall]),split='-'),function(x)x[2])
#c(2,3, 20,21,18,11)]
if(any(baliselistc)==F){
cat('problem chromosome', opt[['chro']], ':', lisc)
cat(unique(Data[,chrohead]))
q(status=2)
}
if(nrow(Data2)==0){
cat("no phenotype ",opt[['pheno']], opt[['file_pheno']]," found in file \n")
q(status=2)
}
Data2Sub<-Data2 #[, c(chrohead,poshead, OrBetaHead,headCi95,pValueHead,phenhead,nvalueHead,freqHead)]

IC<-t(sapply(strsplit(gsub("[", "", sapply(strsplit(as.character(Data2Sub[,headCi95]), split=']',fixed=T),function(x)x[1]),fixed=T), split="-"),function(x){
if(length(x)==2){return(c(as.numeric(x[1]),as.numeric(x[2])))
}else{
return(c(-as.numeric(x[2]),as.numeric(x[3])))
}
}))
IC<-data.frame(lower.cat=IC[,1], upper.cat=IC[,2])
Data2Sub<-cbind(Data2Sub,IC)
Data2Sub$nsample.cat<-sapply(strsplit(as.character(Data2Sub[,nvalueHead]),split="[ ]"),function(x)sum(as.integer(gsub(",", "",grep("[0-9]", x,value=T)))))
Data2Sub$beta.cat<-as.numeric(as.character(Data2Sub[,OrBetaHead]))
balisebeta<-!is.na(Data2Sub$beta.cat) & Data2Sub$beta.cat>=1
Data2Sub$lower.cat[balisebeta]<-log2(Data2Sub$lower.cat[balisebeta])
Data2Sub$upper.cat[balisebeta]<-log2(Data2Sub$upper.cat[balisebeta])

#Data2Sub<-Data2Sub[!is.na(Data2Sub$beta.cat) & !is.na(Data2Sub$lower.cat) & !is.na(Data2Sub$upper.cat) & !is.na(Data2Sub$nsample.cat) ,]

Data2Sub$beta.cat[balisebeta & Data2Sub$beta.cat>=1]<-log2(Data2Sub$beta.cat[balisebeta & Data2Sub$beta.cat>=1])
Data2Sub$risk.allele.af<-as.numeric(as.character(Data2Sub[,freqHead]))
Data2Sub$sd.cat<-(Data2Sub$upper.cat - Data2Sub$beta.cat)/1.96
Data2Sub$sd.cat2<-(Data2Sub$upper.cat - Data2Sub$beta.cat)/1.96*sqrt(Data2Sub$nsample.cat)
Data2Sub$z.cat<-Data2Sub$beta.cat/(Data2Sub$sd.cat)
Data2Sub$h2.cat=computedher(Data2Sub$beta.cat, Data2Sub$sd.cat, Data2Sub$risk.allele.af,Data2Sub$nsample.cat)
Data2Sub$pvalue<-as.numeric(as.character(Data2Sub[,pValueHead]))

write.csv(Data2Sub, file=paste(opt[['out']], '_all.csv',sep=''), row.names=F)

Data2Sub<-Data2Sub[, c(chrohead,poshead,rsHead,riskall,'nsample.cat', 'beta.cat', 'sd.cat', 'z.cat', 'h2.cat', 'pvalue', 'risk.allele.af', 'risk.allele.cat')]
names(Data2Sub)[c(1,2, 3,4)]<-c("chro", "bp", 'rs', "risk_allele")
Data2Sub<-Data2Sub[order(Data2Sub$chro, Data2Sub$bp),]

Data2SubBed<-Data2Sub[,c("chro", "bp", "rs")]
Data2SubBed$bpbefore<-Data2SubBed$bp-1

write.csv(Data2Sub, file=paste(opt[['out']], '_resume.csv',sep=''), row.names=F)
write.table(Data2SubBed[, c("chro", "bpbefore", "bp", "rs")], file=paste(opt[['out']], '.bed',sep=''), row.names=F, col.names=F, sep="\t", quote=F)
write.table(Data2Sub[, c("chro", "bp")], file=paste(opt[['out']], '.pos',sep=''), row.names=F, col.names=F, sep="\t", quote=F)
if(!is.null(opt[['wind']])){
Tmp<-Data2Sub[, c("chro", "bp", "bp","rs")]
Tmp[,2]<-Tmp[,2]-opt[['wind']]*1000
Tmp[,3]<-Tmp[,3]+opt[['wind']]*1000
write.table(Tmp,file=paste(opt[['out']], '_range.bed',sep=''), row.names=F, col.names=F, sep="\t", quote=F)

}
