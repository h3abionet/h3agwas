#!/usr/bin/env Rscript
library(ggplot2)
library(ggplot2)
library("optparse")
getstat<-function(x, head){
  rd<-function(x,roun)round(x, roun)
  if(length(unique(x))==2){
    nmin<-length(which(x==min(x)));ntot<-length(x)
    df<-data.frame(var=head, analyse="n/N (%)",stat=paste(nmin,"/",ntot, " (",round(nmin/ntot*100,1),")",sep=''))
  }else{
   df<-data.frame(var=head, analyse="mean/median (sd, min-max)", stat=paste(rd(mean(x),2),'/', rd(median(x),2), " (",rd(sd(x),2), ", ", rd(min(x),2), '-',rd(max(x),2),")", sep=''))
  }
    return(df)
}

getstatsumm<-function(datapheno, listnewval, filestat){
 Cmt<-1
 for(head in listnewval){
  stat2<-getstat(datapheno[, head], head)
  if(Cmt==1)StatF<-stat2
  else StatF<-rbind(StatF, stat2)
  Cmt<-Cmt+1
 }
 write.csv(StatF, file=filestat, row.names=F, quote=F)
 StatF$var<-as.character(StatF$var);StatF$var[StatF$var=='sex']<-'Female'
 return(StatF)
}


plot_dens<-function(Data, pheno){
mean=mean(Data[,pheno])
fivenum_val=fivenum(Data[,pheno])
resume_stat<-c(mean,fivenum_val)
labels<-paste(c("Mean","0","25","50","75","100"),round(resume_stat,1),sep=':')
print(labels)
ytext<-rep(0.5, length(resume_stat))
ytext[1]<-0.2
data_text<-data.frame(x=resume_stat*1.2, labels=labels, ytext=ytext)
p<-ggplot(Data,aes_string(x=pheno)) +
  geom_density( fill="dodgerblue", alpha=0.5)+
  labs(x=pheno)+
  geom_vline(xintercept=resume_stat, size=1.2, color="red")+
  geom_text(data=data_text,aes(x=x, label=labels, y=ytext))
p
}
checkheader<-function(x, header){
x<-x[!(x %in% header)]
if(length(x)>0){
cat (x, 'not found in header :', paste(header, paste=','))
q('no', 2)
}
}
plot_cmppheno<-function(Data, pheno1, pheno2){
Data1<-data.frame(x=Data[,pheno1],y=Data[,pheno2]);names(Data1)<-c(pheno1, pheno2)
p<-ggplot(Data1, aes_string(x=pheno1, y=pheno2)) + geom_point()
p
}

addresiduals<-function(data, var ,covar, varnewnam,fcttr2=nulltr){
 listvar<-c(var, covar)
 if(length(listvar)>1){
  data2<-data[,c(var, covar)] 
 }else {
 data[,varnewnam]<-fcttr2(data[,var])
 return(data)
 }
 data$num_tmp215<-paste('ind',1:nrow(data),sep='')
 rownames(data2)<-data$num_tmp215
 lmres<-lm(paste(var,"~",  paste(covar, collapse='+')) ,data=data2)
 tmpdf<-data.frame(names(lmres$residuals), fcttr2(lmres$residuals))
 names(tmpdf)<-c('num_tmp215',varnewnam)
 tmpdata<-merge(data,tmpdf,all=T,by="num_tmp215")
 tmpdata
}

invnorm<-function(x)qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))) #inverse norm of residuals



none<-function(x)return(x)

option_list = list(
  make_option(c("--data"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--fam"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--transform_i"), type="character", default='none',
              help="dataset file name", metavar="character"),
  make_option(c("--transform_r"), type="character", default="none",
              help="dataset file name", metavar="character"),
  make_option(c("--covar"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("--residuals"), type="integer", default=1,
              help="dataset file name", metavar="character"),
  make_option(c("--pheno"), type="character", default=1,
              help="dataset file name", metavar="integer"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="output file name [default= %default]", metavar="character")
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

data_all<-read.table(opt[['data']], header=T);headerall<-names(data_all)
datafam<-read.table(opt[['fam']])
data_all<-merge(data_all, datafam, by=c(1,2))

listpheno<-strsplit(opt[['pheno']], split=',')[[1]]
checkheader(listpheno, headerall)
if(!is.null(opt[['covar']])){
listcov<-strsplit(opt[['covar']], split=',')[[1]]
checkheader(listcov, headerall)
}else{
listcov<-c()
}

fcttr1<-get(opt[['transform_i']])
fcttr2<-get(opt[['transform_r']])
out=opt[['out']]
residual<-opt[['residuals']]
if(any(residual %in% c(0,1))==F){
cat ("residual vaslue must be 0 (false) or 1 (true)", residual)
q('no', 2)
}
names(data_all)[c(1,2)]<-c('FID', 'IID')
data_all_tr<-data_all[,c('FID', 'IID', listpheno, listcov)]

newlistpheno<- c()
for(pheno in listpheno){
 newpheno_1<-paste(pheno,'_tr1',sep='')
 data_all_tr[!is.na(data_all_tr[,pheno]),newpheno_1]<-fcttr1(data_all_tr[!is.na(data_all_tr[,pheno]),pheno])
 pp_save<-plot_dens(data_all_tr, pheno)
 ggsave(paste(out,'_',pheno,'_dist.pdf',sep=''),plot=pp_save)
 if(residual>0){
  newpheno<-paste(pheno,'_tr',sep='')
  data_all_tr<-addresiduals(data_all_tr,newpheno_1,listcov, newpheno,fcttr2)
 }else{
  newpheno<-newpheno_1
 }
 psave_cmp<-plot_cmppheno(data_all_tr,pheno, newpheno) 
 ggsave(paste(out,'_',pheno,'_cmp.pdf',sep=''),plot=psave_cmp)
 newlistpheno<-c(newlistpheno,newpheno)
}

if(residual<=0){
writeLines(paste(listcov,collapse=','), con='covar.txt')
}else{
writeLines("", con='covar.txt')
listcov<-c()
}

getstatsumm(data_all_tr, unique(c(newlistpheno, listpheno,listcov)), paste(opt[['out']],"_sumstat.csv",sep=''))

writeLines(paste(newlistpheno,collapse=','), con='pheno.txt')
write.table(data_all_tr[c('FID','IID', newlistpheno, listcov)], row.names=F, col.names=T, sep='\t', quote=F, file=paste(opt[['out']],".pheno",sep=''))
