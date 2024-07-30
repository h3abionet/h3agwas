#!/usr/bin/env Rscript
##  Maintener : jean-Tristan brandenburg
library("optparse")
library(data.table)


option_list = list(
  make_option(c("--vcftools_statref"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--vcftools_stat"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);


read_vcf_tools<-function(File ,head){
data_file<-fread(File,fill=T)
#    CHROM   POS N_ALLELES N_CHR {ALLELE:FREQ}          V6    short_rs polytomie
names(data_file)<-c('CHROM','POS','N_ALLELES','N_CHR',paste('A_FRQ_',1:(ncol(data_file) -4),sep=''))
data_file$short_rs<-paste(data_file$CHROM,data_file$POS,sep=':')
duplicate<-table(data_file$short_rs)
dpnb<-duplicate[duplicate>1]
data_file$polytomie<-(data_file$short_rs %in% names(dpnb)) |data_file$N_ALLELES>2
A1Inf<-strsplit(data_file$A_FRQ_1, split=':')
A2Inf<-strsplit(data_file$A_FRQ_2, split=':')
data_file$A1<-sapply(A1Inf,function(x)x[1])
data_file$A1_frq<-as.numeric(sapply(A1Inf,function(x)x[2]))
data_file$A2<-sapply(A2Inf,function(x)x[1])
data_file$A2_frq<-as.numeric(sapply(A2Inf,function(x)x[2]))
#if(!is.null(head))names(data_file)[3:ncol(data_file)]<-paste(names(data_file)[3:ncol(data_file)],head,sep='_')
return(data_file)
}
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
#opt=list('vcftools_stat'='egt_2017/out_10..vcf.frq','vcftools_statref'='/dataG/cancer/brca/Ark/chr10.frq', out='out')

out=opt[['out']]
data_stat<-read_vcf_tools(opt[['vcftools_stat']],'stat')
data_ref<-read_vcf_tools(opt[['vcftools_statref']],'ref')
data_statc<-data_stat[data_stat$polytomie==F,]
data_refc<-data_ref[data_ref$polytomie==F,]
alldata<-merge(data_statc,data_refc, by=c('CHROM','POS'),all=F, suffixes=c('_stat','_ref'))
#check polytomy
alldata_pol<-merge(data_stat,data_ref, by=c('CHROM','POS'), suffixes=c('_stat','_ref'),all.x=T)
write.table(alldata_pol, file=paste(out,'.pol.tsv',sep=''))
alldata_pol<-na.omit(alldata_pol)
alldata_pol2<-alldata_pol[alldata_pol$polytomie_stat | alldata_pol$polytomie_ref,]
write.csv(alldata_pol2, file=paste(out,'.pol.csv',sep=''))
uniq_dupl<-unique(alldata_pol2[,c('CHROM','POS')])
balise<-((alldata$A1_stat==alldata$A1_ref) &(alldata$A2_stat==alldata$A2_ref)) |  ((alldata$A1_stat==alldata$A2_ref) &(alldata$A2_stat==alldata$A1_ref))
write.table(alldata, file=paste(out,'_mergebeforefilt.csv',sep=''))
alldatafilter=alldata[balise,]
# balise 
balrev<-alldatafilter$A1_stat !=alldatafilter$A1_ref
alldatafilter$A1_frq_stat[balrev]<- 1 - alldatafilter$A1_frq_stat[balrev]
alldatafilter$diff_freq<-abs(alldatafilter$A1_frq_ref-alldatafilter$A1_frq_stat)

#    CHROM   POS N_ALLELES_stat N_CHR_stat A_FRQ_1_stat A_FRQ_2_stat
#      short_rs_stat polytomie_stat A1_stat A1_frq_stat A2_stat A2_frq_stat
#      N_ALLELES_ref N_CHR_ref A_FRQ_1_ref A_FRQ_2_ref short_rs_ref polytomie_ref
#      A1_ref A1_frq_ref A2_ref A2_frq_ref
alldatafilter2<-alldatafilter[,c('CHROM','POS','N_CHR_stat','N_CHR_ref','A1_ref','A2_ref', 'A1_frq_ref', 'A2_frq_ref')]
write.table(alldatafilter2, file=paste(out, '_clean.tsv',sep=''),quote=F, sep='\t',col.names=T, row.names=F)
# CHROM,POS
#resumestat<-data.frame(nbpi_stat=nrow(unique(data_stat[,c('CHROM','POS')])), nbpi_ref=nrow(unique(data_ref[,c('CHROM','POS')])),nbbpclean_stat=nrow(data_statc),nbbpclean_ref=nrow(data_refc), common_pos=nrow(alldata), poly_dup_nb=nrow(uniq_dupl), nb_allele_bad=length(which(!balise)), r2freq=cor(alldatafilter2$A1_frq_stat,alldatafilter2$A1_frq_ref))
#write.csv(resumestat, quote=F, row.names=F, file=paste(out, '_resume.csv',sep=''))

