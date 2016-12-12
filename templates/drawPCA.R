#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
eigval     <- read.table("$EIGVALS")
eigvec1.pc <- round(eigval[1,1]/sum(eigval)*100,digits=2)
eigvec2.pc <- round(eigval[2,1]/sum(eigval)*100,digits=2)
eigvec     <- read.table("$EIGVECS")
pdf("$OUTPUT")
par(xpd = T, mar = par()\$mar + c(0,0,0,7))
plot(eigvec[,3], eigvec[,4], 
   xlab=paste("PCA 1\n",eigvec1.pc, "% of observed variation", sep=""),
   ylab=paste("PCA 2\n",eigvec2.pc, "% of observed variation", sep=""),
   col=c("blue","hotpink"),cex=1,cex.lab=1.1,cex.axis=1.1)

legend(0.05,0, pch=c(1,15), col=c("blue","hotpink"), legend = c("Case","Control"), cex=1.5)
par(mar=c(5, 4, 4, 2) + 0.1)