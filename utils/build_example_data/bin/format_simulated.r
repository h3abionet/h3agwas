#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-g", "--gc_her"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-b", "--bfile"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-y", "--clump_p1"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-z", "--clump_p2"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-r", "--clump_r2"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-k", "--clump_kb"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-i", "--nb_snp"), type="integer", default=-1,
              help="dataset file name", metavar="character"),
  make_option(c("--used_effect"), type="integer", default=1,
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

Datainfo<-read.csv(opt[['gc_her']])

bimfile<-read.table(paste(opt[['bfile']],'.bim',sep=''), sep='\t')
names(bimfile)<-c('chro', 'rs_bim', 'cm', 'bp','alt', 'ref')
print(head(bimfile))
print(head(Datainfo))
bimfile$chr<-as.character(bimfile$chr)
allinfo<-merge(Datainfo, bimfile, by=c('chro', 'bp'))
if(nrow(allinfo)==0){
cat('not commomn position between position rs in gwas catalog and genotype\n')
q(status=2)
}

allinfo$z.cat.new<-allinfo$z.cat
allinfo$riskall<-sapply(strsplit(as.character(allinfo$risk_allele),split='-'),function(x)x[2])
balise<-allinfo$riskall==allinfo$ref
allinfo$z.cat.new[balise]<- -allinfo$z.cat.new[balise]
allinfo$beta.cat[balise]<- -allinfo$beta.cat[balise]
allinfo$A2<-allinfo$ref
allinfo$A2[balise]<-allinfo$alt[balise]
allinfo$A1<-allinfo$alt
allinfo$A1[balise]<-allinfo$ref[balise]
## clump data
#CHR	Chromosome code. Requires --map.
#SNP	Variant identifier.
#BP	Base-pair coordinate. Requires --map.
#A1	Allele 1 (usually minor)
#A2	Allele 2 (usually major)
#FRQ	Allele 1 frequency
#INFO	R-squared quality metric/information content
#'BETA'/'OR'	Regression coefficient (for quantitative traits) or odds ratio
#SE	Standard error of effect (not odds ratio) estimate
#P	Association test p-value
#[c('chro','rs_bim', 'bp','A1', 'A2', 'risk.allele.af', 'beta.cat', 'se.cat', 'pvalue')]
#names(allinfo)
#print(c('chro','rs_bim', 'bp','A1', 'A2', 'risk.allele.af', 'beta.cat', 'se.cat', 'pvalue')[!(c('chro','rs_bim', 'bp','A1', 'A2', 'risk.allele.af', 'beta.cat', 'se.cat', 'pvalue') %in% names(allinfo))])
dataformatplk<-allinfo[,c('chro','rs_bim', 'bp','A1', 'A2', 'risk.allele.af', 'beta.cat', 'sd.cat', 'pvalue')]
names(dataformatplk)<-c("CHR", "SNP", "BP","A1", "A2", "FRQ", "BETA", "SE", "P")
write.table(dataformatplk,file=paste(opt[['out']], ".assoc",sep=''),sep="\t", row.names=F,col.names=T,quote=F) 
system(paste("plink -bfile ",opt[['bfile']]," --clump ", opt[['out']], ".assoc", " --clump-p1 ", opt[['clump_p1']], " --clump-p2 ", opt[['clump_p2']]," --clump-r2 ", opt[['clump_r2']], " --clump-kb ", opt[['clump_kb']], ' -out ', opt[['out']], '_clump' ,sep=''))
Dataclump<-read.table(paste(opt[['out']], '_clump.clumped',sep=''), header=T)

allinfo<-allinfo[allinfo$rs_bim %in% Dataclump$SNP ,]

rsbest<-aggregate(z.cat.new~rs_bim, allinfo,function(x)x[which.min(abs(x))])
#rsbest$z.cat.new<-rsbest$z.cat.new*opt[['increase']]
nb_snp=opt[['nb_snp']]
if(nb_snp>0 & nb_snp<=nrow(rsbest))rsbest<-rsbest[sample(1:nrow(rsbest), nb_snp),]

if(opt[['used_effect']]==0)rsbest[,2]<-rnorm(nrow(rsbest))
names(rsbest)<-c('rsid', 'z.simulated')
allinfo<-merge(allinfo,rsbest,by.x='rs_bim', by.y='rsid',all=F)
write.table(allinfo, sep='\t', quote=F, file=paste(opt[['out']],'.infors',sep=''), row.names=F, col.names=T)
write.table(rsbest, sep='\t', quote=F, file=opt[['out']], row.names=F, col.names=F)


