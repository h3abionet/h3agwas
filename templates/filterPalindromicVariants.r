#!/usr/bin/env Rscript

#  FILTER PALINDROMIC VARIANTS
#  ===========================
#
#  This script marks and filters out palindromic [A/T]/[C/G]
#  variants from the PLINK .bim file. These variants will be
#  excluded prior to phasing
#
#
##############################################################

bimFile = read.table("${params.outputDir}${bim}", header=F)

markedPalindromes = which((bimFile\$V5=="A" & bimFile\$V6=="T") | +
			  (bimFile\$V5=="T" & bimFile\$V6=="A") | +
			  (bimFile\$V5=="G" & bimFile\$V6=="C") | +
			  (bimFile\$V5=="C" & bimFile\$V6=="G"))

selectedPalindromes = bimFile[markedPalindromes,]

write.table(selectedPalindromes\$V2, 
	    file="${params.cohortName}.palindromic.snvs.txt", 
	    row.names=F, 
	    col.names=F, 
	    quote=F)
