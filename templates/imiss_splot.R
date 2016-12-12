#!/usr/bin/env Rscript
#Load SNP frequency file and generate histogram
b.frq <- read.table("$input",header=T)
pdf("$output")
plot(ecdf(b.frq\$F_MISS),xlim=c(0,0.10),ylim=c(0,1),pch=20, main="Individual Missingness Distribution", xlab="Missingness Frequency", ylab="Fraction of SNPs",col="blue",axes=T)
