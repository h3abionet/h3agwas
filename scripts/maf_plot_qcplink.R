#!/usr/bin/env Rscript
#Load SNP frequency file and generate cumulative freequency distribution
args <- commandArgs(TRUE)
b.frq <- read.table(args[1],header=T)
pdf(args[2])
plot(ecdf(b.frq$MAF), xlim=c(0,0.10),ylim=c(0,1),pch=20, main="MAF cumulative distribution",xlab="Minor allele frequency (MAF)", ylab="Fraction of SNPs",axes=T)
