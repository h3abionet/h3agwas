#!/usr/bin/env Rscript
#Load SNP frequency file and generate histogram
args <- commandArgs(TRUE)
b.frq <- read.table(args[1],header=T)
pdf(args[2])
plot(ecdf(b.frq$F_MISS),xlim=c(0,0.10),ylim=c(0,1),pch=20, main="SNP Missingness Distribution", xlab="Missingness Frequency", ylab="Fraction of SNPs",col="blue",axes=T)
