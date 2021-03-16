#!/usr/bin/env Rscript
library(data.table)
args = commandArgs(trailingOnly=TRUE)
databim=fread(args[1],header=F)
databim$chrpos<-paste(databim$V1, databim$V4)
tmp=table(databim$chrpos)
tmp2<-table(as.character(databim$V2))
balise<-(as.character(databim$V5) %in% c("A", "T", "C", "G")) & (as.character(databim$V6) %in% c("A", "T", "C", "G"))  & (databim$chrpos %in% names(tmp[tmp==1])) & (databim$V2 %in% names(tmp2[tmp2==1]))
listrstodel=as.character(databim$V2[!balise])
#system("plink --keep-allele-order --make-bed --bfile $plk --out $out -maf 0.0000000000000000001 --exclude rstodel --threads ${params.max_plink_cores}
if(length(listrstodel)>0){
writeLines(as.character(databim$V2[!balise]), con=args[2])
}else{
writeLines("nors", con=args[2])
}

