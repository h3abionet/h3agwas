#!/usr/bin/Rscript
stat_box_data <- function(x) {
  upper_limit=min(x, na.rm=T) - .05*min(x, na.rm=T)
  return( 
    data.frame(
      y = 0.95 * upper_limit,
      label = paste('count =', 
                    format(length(x), big.mark = ",", decimal.mark = ".", scientific = FALSE), 
                    '\n',
                    'mean =', 
                    format(round(mean(x), 3), big.mark = ",", decimal.mark = ".", scientific = FALSE))
    )
  )
}
ListPackNeed=c("ggpubr", "optparse", "gridExtra")
ListPackNeedIn<-ListPackNeed[!ListPackNeed %in% rownames(installed.packages())]
for(pack in ListPackNeedIn)if(pack %in% rownames(installed.packages()) == FALSE) {install.packages(pack, lib= Sys.getenv("R_LIBS_USER"), repos='http://cran.us.r-project.org')}
library("ggpubr")
library("optparse")
require(gridExtra)


#!/usr/bin/Rscript
as.characterspe<-function(x,round=2){
y<-rep(NA, length(x))
bal=!is.na(x) & x>10**-round
y[bal]<-as.character(format(round(x[bal],round), scientific=F))
y[!bal]<-as.character(format(x[!bal], scientific=T,digits=round))
y
}

GetGPlotBox<-function(datamerg, pheno, title){
UnRes<-unique(datamerg$GENO)
if(length(UnRes)==3){
my_comparisons <- list(c(UnRes[1], UnRes[3]), c(UnRes[2], UnRes[3]), c(UnRes[1], UnRes[2]) )
reskr<-kruskal.test(as.formula(paste(pheno,"~ GENO")), data = datamerg)
p <- ggboxplot(datamerg, x = "GENO", y = pheno, color = "GENO", palette = "jco", title=title,subtitle=paste("Kruskal-Wallis test, p = ", as.characterspe(reskr$p.value),sep=""))+stat_compare_means(comparisons=my_comparisons) + stat_summary(
    fun.data = stat_box_data,
    geom = "text")
}else{
my_comparisons <- list(c(UnRes[1], UnRes[2]))
p <- ggboxplot(datamerg, x = "GENO", y = pheno, color = "GENO", palette = "jco") + stat_compare_means(comparisons=my_comparisons) + stat_summary(
    fun.data = stat_box_data, 
    geom = "text"
  )
}
p
}
getresidual<-function(datamerge, pheno, cov){
datamerg$FIDIID<-paste(datamerg[,1], datamerg[,2])
rownames(datamerg)<-datamerg$FIDIID
Res<-glm(as.formula(paste(pheno,"~", paste(cov, collapse="+"))), data=datamerg)$residuals
datamergres=merge(datamerg, data.frame(FIDIID=names(Res), residuals=Res), by="FIDIID")
return(datamergres)
}


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
  make_option(c("-g", "--gxe"), type="character",
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
gxe<-c()
if(!is.null(opt[['gxe']]))gxe=opt[['gxe']]
if(!is.null(opt[['cov']]))cov=gsub(" ","",strsplit(opt[['cov']],split=",")[[1]]) 

datamerg<-merge(data_pheno[,c(names(data_pheno)[1:2],pheno,cov, gxe)], data_ped)
datamerg<-na.omit(datamerg)

if(length(cov)>0)datamergres<-getresidual(datamerge, pheno, cov)

if(length(gxe)==0){
p<-GetGPlotBox(datamerg,pheno, pheno)
if(length(cov)>0){
p2<-GetGPlotBox(datamergres,pheno, paste("res ",paste(pheno,"~", paste(cov, collapse="+"))))
gg<-arrangeGrob(p, p2, ncol=2)
}else{
gg<-p
}
}else{
stat1=unique(datamerg[,gxe])[1]
stat2=unique(datamerg[,gxe])[2]
p<-GetGPlotBox(datamerg[datamerg[, gxe]==stat1,],pheno, paste(pheno,", ", gxe, ": ", stat1,sep=""))
p2<-GetGPlotBox(datamerg[datamerg[, gxe]==stat2,],pheno, paste(pheno,", ", gxe, ": ", stat2,sep=""))
if(length(cov)>0){
p3<-GetGPlotBox(datamergres[datamergres[, gxe]==stat1,],pheno, paste("res ",pheno,", ", gxe, ": ", stat1,sep=""))
p4<-GetGPlotBox(datamergres[datamergres[, gxe]==stat2,],pheno, paste("res ",pheno,", ", gxe, ": ", stat2,sep=""))
gg<-arrangeGrob(p, p2, p3,p4, ncol=2, nrow=2)
}else{
gg<-arrangeGrob(p, p2, ncol=2)
}
}
ggsave(opt[['out']],gg,width=7, height=7)



