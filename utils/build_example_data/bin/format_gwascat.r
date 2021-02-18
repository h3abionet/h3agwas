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
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-p", "--pheno"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-c", "--chro"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#bin 	590	smallint(5) unsigned 	range 	Indexing field to speed chromosome range queries.
#chrom 	chr1	varchar(255) 	values 	Reference sequence chromosome or scaffold
#chromStart 	768252	int(10) unsigned 	range 	Start position in chromosome
#chromEnd 	768253	int(10) unsigned 	range 	End position in chromosome
#name 	rs2977608	varchar(255) 	values 	ID of SNP associated with trait
#pubMedID 	31969693	int(10) unsigned 	range 	PubMed ID of publication of the study
#author 	Coleman JRI	varchar(255) 	values 	First author of publication
#pubDate 	2020-01-23	varchar(255) 	values 	Date of publication
#journal 	Mol Psychiatry	varchar(255) 	values 	Journal of publication
#title 	Genome-wide gene-environmen...	varchar(1024) 	values 	Title of publication
#trait 	Major depressive disorder	varchar(255) 	values 	Disease or trait assessed in study
#initSample 	29,475 European ancestry ca...	longblob 	  	Initial sample size
#replSample 	NA	longblob 	  	Replication sample size
#region 	1p36.33	varchar(255) 	values 	Chromosome band / region of SNP
#genes 	NR	longblob 	  	Reported Gene(s)
#riskAllele 	rs2977608-A	longblob 	  	Strongest SNP-Risk Allele
#riskAlFreq 	0.259029	varchar(255) 	values 	Risk Allele Frequency
#pValue 	8E-6	varchar(255) 	values 	p-Value
#pValueDesc 	 	varchar(255) 	values 	p-Value Description
#orOrBeta 	1.0729614	varchar(255) 	values 	Odds ratio or beta
#ci95 	[1.04-1.1]	varchar(255) 	values 	95% Confidence Interval
#platform 	Affymetrix [7791636] (imputed)	varchar(255) 	values 	Platform and [SNPs passing QC]
#cnv 	N	enum('Y', 'N') 	values 	Y if Copy Number Variant

lisp<-unlist(strsplit(opt[['pheno']],','))
lisc<-unlist(strsplit(opt[['chro']],','))
Data<-read.csv(opt[['file']], header=F, sep='\t')

#lisp="Type 2 diabetes";lisc=c("22","21")
format="USCS"
if(format=="USCS"){
chrohead="chrom";phenhead="trait";poshead="chromEnd";OrBetaHead="orOrBeta";headCi95="ci95";pValueHead="pValue";nvalueHead="initSample";freqHead="riskAlFreq";rsHead="name";riskall="riskAllele"
names(Data)<-c("bin", "chrom", "chromStart", "chromEnd", "name", "pubMedID", "author", "pubDate", "journal", "title", "trait", "initSample","replSample","region", "genes", "riskAllele", "riskAlFreq", "pValue", "pValueDesc", "orOrBeta", "ci95","platform", "cnv")
}
Data[,chrohead]<-gsub('chr', '',as.character(Data[,chrohead]))
Data2<-Data[tolower(as.character(Data[,phenhead])) %in% tolower(lisp) & Data[,chrohead] %in% lisc,]
#c(2,3, 20,21,18,11)]
if(nrow(Data2)==0){
cat("no phenotype ",opt[['pheno']]," found in file \n")
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
balisebeta<-!is.na(Data2$beta.cat) & Data2Sub$beta.cat>=1
Data2Sub$lower.cat[balisebeta]<-log2(Data2Sub$lower.cat[balisebeta])
Data2Sub$upper.cat[balisebeta]<-log2(Data2Sub$upper.cat[balisebeta])
Data2Sub<-Data2Sub[!is.na(Data2Sub$beta.cat) & !is.na(Data2Sub$lower.cat) & !is.na(Data2Sub$upper.cat) & !is.na(Data2Sub$nsample.cat) ,]

Data2Sub$beta.cat[Data2Sub$beta.cat>=1]<-log2(Data2Sub$beta.cat[Data2Sub$beta.cat>=1])
Data2Sub$risk.allele.af<-as.numeric(as.character(Data2Sub[,freqHead]))
Data2Sub$sd.cat<-(Data2Sub$upper.cat - Data2Sub$beta.cat)/1.96
Data2Sub$sd.cat2<-(Data2Sub$upper.cat - Data2Sub$beta.cat)/1.96*sqrt(Data2Sub$nsample.cat)
Data2Sub$z.cat<-Data2Sub$beta.cat/(Data2Sub$sd.cat)
Data2Sub$h2.cat=computedher(Data2Sub$beta.cat, Data2Sub$sd.cat, Data2Sub$risk.allele.af,Data2Sub$nsample.cat)
Data2Sub$pvalue<-as.numeric(as.character(Data2Sub[,pValueHead]))
Data2Sub<-Data2Sub[!is.na(Data2Sub$h2.cat) & !is.na(Data2Sub[,pValueHead]),]
Data2Sub<-Data2Sub[order(Data2Sub$h2.cat),]
Data2Sub$order<-1:nrow(Data2Sub)
Good<-aggregate(as.formula(paste("order~",chrohead,"+",poshead,sep="")), data=Data2Sub, min)$order
Data2Sub<-Data2Sub[Data2Sub$order %in% Good,]

write.csv(Data2Sub, file=paste(opt[['out']], '_all.csv',sep=''), row.names=F)

Data2Sub<-Data2Sub[, c(chrohead,poshead,rsHead,riskall,'nsample.cat', 'beta.cat', 'sd.cat', 'z.cat', 'h2.cat', 'pvalue', 'risk.allele.af')]
names(Data2Sub)[c(1,2, 3,4)]<-c("chro", "bp", 'rs', "risk_allele")
Data2Sub<-Data2Sub[order(Data2Sub$chro, Data2Sub$bp),]
write.csv(Data2Sub, file=paste(opt[['out']], '_resume.csv',sep=''), row.names=F)

Data2SubBed<-Data2Sub[,c("chro", "bp", "rs")]
Data2SubBed$bpbefore<-Data2SubBed$bp-1

write.table(Data2SubBed[, c("chro", "bpbefore", "bp", "rs")], file=paste(opt[['out']], '.bed',sep=''), row.names=F, col.names=F, sep="\t", quote=F)
write.table(Data2Sub[, c("chro", "bp")], file=paste(opt[['out']], '.pos',sep=''), row.names=F, col.names=F, sep="\t", quote=F)
