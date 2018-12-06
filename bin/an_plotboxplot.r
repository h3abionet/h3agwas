#!/usr/bin/Rscript

#if("ggpubr" %in% rownames(installed.packages()) == FALSE) {install.packages("ggpubr", lib=.libPaths()[2])}
#!/usr/bin/Rscript
as.characterspe<-function(x,round=2){
y<-rep(NA, length(x))
bal=!is.na(x) & x>10**-round
y[bal]<-as.character(format(round(x[bal],round), scientific=F))
y[!bal]<-as.character(format(x[!bal], scientific=T,digits=round))
y
}

library("ggpubr")
library("optparse")
require(gridExtra)


option_list = list(
  make_option(c("-d", "--data"), type="character",
              help="data files with phenotype and FID", metavar="character"),
  make_option(c("-p", "--ped"), type="character", 
              help="ped file contains genotype for each ", metavar="character"),
  make_option(c("-e", "--pheno"), type="character", 
              help="ped file contains genotype for each ", metavar="character"),
  make_option(c("-c", "--cov"), type="character", 
              help="ped file contains genotype for each ", metavar="character"),
  make_option(c("-t", "--type_out"), type="character", default="pdf",
              help="type output file name [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out",
              help="output file name [default= %default]", metavar="character")
);

args = commandArgs(trailingOnly=TRUE)
if(length(args)<2){
opt=list(ped="test/plk_rs190712692.ped",data="test/All.pheno",out="test_rs1.pdf", pheno="cholesterol_1_qc", cov="sex,age")
}else{
opt = OptionParser(option_list=option_list);
opt = parse_args(opt);
}
data_ped=read.table(opt[["ped"]], sep="\t")
data_pheno=read.table(opt[["data"]], header=T)
data_ped<-data_ped[, c(1,2,ncol(data_ped))]
names(data_ped)<-c("FID","IID","GENO")
pheno<-opt[['pheno']]
cov<-c()
if(!is.null(opt[['cov']]))cov=gsub(" ","",strsplit(opt[['cov']],split=",")[[1]]) 
datamerg<-merge(data_pheno[,c(names(data_pheno)[1:2],pheno,cov)], data_ped)
datamerg<-na.omit(datamerg)
UnRes<-unique(datamerg$GENO)
if(length(UnRes)==3){
my_comparisons <- list(c(UnRes[1], UnRes[3]), c(UnRes[2], UnRes[3]), c(UnRes[1], UnRes[2]) )
#datamerg$GENO<-as.character(datamerg$GENO)
reskr<-kruskal.test(as.formula(paste(pheno,"~ GENO")), data = datamerg)
p <- ggboxplot(datamerg, x = "GENO", y = pheno, color = "GENO", palette = "jco", title=pheno,subtitle=paste("Kruskal-Wallis test, p = ", as.characterspe(reskr$p.value),sep="")) + stat_compare_means(comparisons=my_comparisons) 
}else{
my_comparisons <- list(c(UnRes[1], UnRes[2]))
p <- ggboxplot(datamerg, x = "GENO", y = pheno, color = "GENO", palette = "jco") + stat_compare_means(comparisons=my_comparisons)
}
#  Add p-value
if(length(cov)>0){
datamerg$FIDIID<-paste(datamerg[,1], datamerg[,2])
rownames(datamerg)<-datamerg$FIDIID
Res<-glm(as.formula(paste(pheno,"~", paste(cov, collapse="+"))), data=datamerg)$residuals
datamergres=merge(datamerg, data.frame(FIDIID=names(Res), residuals=Res), by="FIDIID")

if(length(UnRes)==3){
my_comparisons <- list(c(UnRes[1], UnRes[3]), c(UnRes[2], UnRes[3]), c(UnRes[1], UnRes[2]) )
reskr<-kruskal.test(as.formula(paste(pheno,"~ residuals")), data = datamergres)
#datamerg$GENO<-as.character(datamerg$GENO)
p2 <- ggboxplot(datamergres, x = "GENO", y = "residuals", color = "GENO", palette = "jco", title=paste(paste(pheno,"~", paste(cov, collapse="+"))) ,subtitle=paste("Kruskal-Wallis test, p = ", as.characterspe(reskr$p.value))) + stat_compare_means(comparisons=my_comparisons)
}else{
my_comparisons <- list(UnRes)
p2 <- ggboxplot(datamergres, x = "GENO", y = "residuals", color = "GENO", palette = "jco", title=paste("res ",paste(pheno,"~", paste(cov, collapse="+")))) + stat_compare_means(comparisons=my_comparisons)
}
#pdf(opt[['out']],width=7, height=7)
grid.arrange(p, p2, ncol=2)
ggsave(opt[['out']],width=7, height=7)
#dev.off()
}else{

#pdf(opt[['out']])
p
ggsave(opt[['out']],width=7, height=7)

#dev.off()
}


