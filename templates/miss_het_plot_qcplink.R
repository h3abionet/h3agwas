#!/usr/bin/env Rscript
#--INSPECT MISSINGNESS PATTERNS--#

#IMPORT PLINK FILES WITH MISSINGNESS INFORMATION
#requires the files qcplink_miss.imiss and qcplink_het.het to be present in the script folder

args <- commandArgs(TRUE)
imiss <- read.table("$imiss",header=T)
het <- read.table("$het",header=T)

#CHECK THAT THE	PROPORTION OF MISSING GENOTYPES	IS NOT O
#NOTE: IF F_MISS IS ZERO THEN WE ONLY PLOT MEAN HETEROZYGOSITY

if (!(min(imiss\$F_MISS) == 0 && max(imiss\$F_MISS) == 0)) {


   #CALCULATE CALL RATE, LOG10(F_FMISS) and mean heterozygosity
   imiss\$CALL_RATE <- 1-imiss\$F_MISS
   imiss\$logF_MISS = log10(imiss[,6])
   het\$meanHet = (het\$N.NM. - het\$O.HOM.)/het\$N.NM.
   het\$meanHet <- ifelse(het\$meanHet=="NaN", c(0),c(het\$meanHet))
   imiss.het <- merge(het,imiss,by=c("FID","IID"))

   #Print Heterozygosity cutoffs
   print(paste("cut_het_low: heterozygosity_mean - 3sd is ",sprintf("%.3f", mean(het\$meanHet)-(3*sd(het\$meanHet)), sep="")))
   print(paste("cut_het_high: heterozygosity_mean + 3sd is ",sprintf("%.3f", mean(het\$meanHet)+(3*sd(het\$meanHet))), sep=""))

   #GENERATE CALL RATE BY HETEROZYGOSITY PLOT
   colors  <- densCols(imiss\$logF_MISS,het\$meanHet)
   pdf("$pairs")
   #plot(imiss\$logF_MISS,het\$meanHet, col=colors, xlim=c(-3,0),ylim=c(0.26,0.35),pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate", axes=F)
   plot(imiss\$logF_MISS,het\$meanHet, col=colors, xlim=c(-3,0),ylim=c(0,0.5), pch=20, xlab="Proportion of missing genotypes", ylab="Heterozygosity rate", axes=F)
   #axis(2,at=c(0.26,0.27,0.28, 0.29,0.3,0.31,0.32,0.33,0.34,0.35),tick=T)
   axis(2,at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5),tick=T)
   axis(1,at=c(-3,-2,-1,0),labels=c(0.001,0.01,0.1,1))
   #Heterozygosity thresholds (Horizontal Line)
   abline(h=mean(het\$meanHet)-(3*sd(het\$meanHet)),col="RED",lty=2)
   abline(h=mean(het\$meanHet)+(3*sd(het\$meanHet)),col="RED",lty=2)
   #Missing Data Thresholds (Vertical Line)
   abline(v=-1.30103, col="BLUE", lty=2) #THRESHOLD=0.07
   abline(v=-1.522879, col="RED", lty=2) #THRESHOLD=0.05
   
} else {
  
    het\$meanHet = (het\$N.NM. - het\$O.HOM.)/het\$N.NM.
    het\$meanHet <- ifelse(het\$meanHet=="NaN", c(0),c(het\$meanHet))

    #Print Heterozygosity cutoffs
    print(paste("cut_het_low: heterozygosity_mean - 3sd is ",sprintf("%.3f", mean(het\$meanHet)-(3*sd(het\$meanHet)), sep="")))
    print(paste("cut_het_high: heterozygosity_mean + 3sd is ",sprintf("%.3f", mean(het\$meanHet)+(3*sd(het\$meanHet))), sep=""))
    
    pdf("$meanhet")	
    plot(het\$meanHet)
    abline(h=mean(het\$meanHet)-(3*sd(het\$meanHet)),col="RED",lty=2)
    abline(h=mean(het\$meanHet)+(3*sd(het\$meanHet)),col="RED",lty=2)}
