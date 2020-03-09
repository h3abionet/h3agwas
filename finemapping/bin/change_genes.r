#!/usr/bin/env Rscript

library(data.table)
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
filegenes="/spaces/jeantristan/GWAS/Ressource/gencode.v19.annotation.gtf"
}else{
filegenes=args[1]
}
DataGenes<-as.data.frame(fread(filegenes, sep="\t", header=F))
DataGenes<-DataGenes[DataGenes$V3 %in% c('CDS', 'gene'),]
DataGenes$PosGene<-1:nrow(DataGenes)
DataGenes$GeneName<-sapply(strsplit(DataGenes$V9,split=";"),function(x)gsub("\"","",gsub(" ","",gsub("gene_name","",grep("gene_name",x, value=T)))))
DataGenes2<-DataGenes[,c(1,3,4,5,11)]
names(DataGenes2)<-c("CHR", "Type","BEGIN", "END","GENE")
DataGenes2$CHR<-gsub('chr','',DataGenes2$CHR)
write.table(DataGenes2, sep="\t", row.names=F, file="gencode.v19.genes", col.names=T)



