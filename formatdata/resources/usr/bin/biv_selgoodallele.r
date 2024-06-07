#!/usr/bin/env Rscript
library(data.table)
args = commandArgs(trailingOnly=TRUE)

bim=args[1]

justatgc=T
if(is.na(args[3])){
	justatgc=T
}else{
if(args[3]=='1'){
justatgc=F
}
}


databim=fread(bim,header=F)
databim$chrpos<-paste(databim$V1, databim$V4)
tmp=table(databim$chrpos)
tmp2<-table(as.character(databim$V2))
if(justatgc){
balise<-(as.character(databim$V5) %in% c("A", "T", "C", "G")) & (as.character(databim$V6) %in% c("A", "T", "C", "G"))  & (databim$chrpos %in% names(tmp[tmp==1])) & (databim$V2 %in% names(tmp2[tmp2==1]))
}else{
balise<-(databim$chrpos %in% names(tmp[tmp==1])) & (databim$V2 %in% names(tmp2[tmp2==1]))
}
listrstodel=as.character(databim$V2[!balise])
#system("plink --keep-allele-order --make-bed --bfile $plk --out $out -maf 0.0000000000000000001 --exclude rstodel --threads ${params.max_plink_cores}
if(length(listrstodel)>0){
writeLines(listrstodel, con=args[2])
}else{
writeLines("nors", con=args[2])
}

