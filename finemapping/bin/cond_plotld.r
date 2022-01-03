#!/usr/bin/env Rscript
#library(gaston)
library("optparse")
snp.sq <- function(i, j, ld, write.ld.val, color, polygon.par, cex.ld) {
  if(i == j) return(); # pas de plot d'un SNP avec lui même
  if(j < i) { tmp <- j; j <- i; i <- tmp } # i < j
  d <- (j-i)/2
  cx <- i+d
  cy <- -d
  do.call(polygon, c(list(x = cx+c(-1,0,1,0)/2, y = cy + c(0,1,0,-1)/2, col = color), polygon.par))
  if(write.ld.val) text(cx, cy, ld, cex = cex.ld)
}



LD.plot2<-function (LD, snp.positions, max.dist = Inf, depth = nrow(LD), 
    graphical.par = list(mar = c(0, 0, 0, 0)), cex.ld, cex.snp, 
    polygon.par = list(border = "white"), color.scheme = function(ld) rgb(1, 
        1 - abs(ld), 1 - abs(ld)), write.snp.id = TRUE, write.ld = function(ld) sprintf("%.2f", 
        ld), draw.chr = TRUE, above.space = 1 + 2 * write.snp.id + 
        draw.chr, below.space = 1, pdf.file, finalize.pdf = TRUE, write.snp.id.col=list()) 
{
    n <- nrow(LD)
    positions <- if (missing(snp.positions)) 
        rep(0, n)
    else snp.positions
    graph.depth <- 0
    for (i in seq(1, n - 1)) for (j in seq(i + 1, min(n, i + 
        depth))) {
        if (positions[j] - positions[i] > max.dist) 
            next
        graph.depth <- max(graph.depth, (j - i)/2)
    }
    if (!missing(pdf.file)) 
        pdf(pdf.file, width = n/2, height = (graph.depth + above.space + 
            below.space)/2)
    do.call(par, graphical.par)
    write.ld.val <- ifelse(is.null(write.ld), FALSE, TRUE)
    plot(0, 0, xlim = c(0.5, n + 0.5), ylim = c(-graph.depth - 
        below.space, above.space), type = "n", asp = 1, xaxt = "n", 
        yaxt = "n", bty = "n", xlab = "", ylab = "")
    ld <- ""
    if (missing(cex.ld) & write.ld.val) 
        cex.ld = 0.7 * par("cex")/max(strwidth(sapply(LD, write.ld)))
    for (i in seq(1, n - 1)) for (j in seq(i + 1, min(n, i + 
        depth))) {
        if (positions[j] - positions[i] > max.dist) 
            next
        if (write.ld.val) 
            ld <- write.ld(LD[i, j])
        snp.sq(i, j, ld, write.ld.val, color.scheme(LD[i, j]), 
            polygon.par, cex.ld)
    }
    if (write.snp.id) {
        rs.h <- strheight(rownames(LD))
        if (missing(cex.snp)) 
            cex.snp <- 0.25 * par("cex")/max(rs.h)
        col=rep('black', n)
        rsname=rownames(LD)
        if(length(write.snp.id.col)>0){
          for(col_cmp in names(write.snp.id.col)){
           col[rsname %in% write.snp.id.col[[col_cmp]]]<-col_cmp 
          }
        }
        text(1:n, 0, rownames(LD), srt = 90, cex = cex.snp, adj = c(0, 
            0.5), col=col)
    }
    if (!missing(snp.positions) & draw.chr) {
        if (write.snp.id) 
            a <- max(strwidth(rownames(LD), cex = cex.snp))
        else a <- 0
        pos <- 1.5 + (n - 2) * (snp.positions - snp.positions[1])/diff(range(snp.positions))
        segments(1:n, a + 0.25, pos, a + 1.5)
        rect(1.5, a + 1.5, n - 0.5, a + 1.75)
        segments(pos, a + 1.5, pos, a + 1.75)
    }
    if (!missing(pdf.file) & finalize.pdf) 
        dev.off()
}

 
option_list = list(
  make_option(c("--ld"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--bim"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--col_rs"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--pos_ref"), type="integer", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("--out"), type="character", default="out.svg", 
              help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

ld=opt[['ld']]
bim=opt[['bim']]
out=opt[['out']]


m <- as.matrix(read.csv(ld, sep="\t", header=F))
DataBim<-read.table(bim)
colnames(m)<-DataBim$V2
rownames(m)<-DataBim$V2
write.snp.id.col=list()
if(!is.null(opt[['col_rs']])){
DataCol<-read.table(opt[['col_rs']])
DataCol2<-merge(DataCol, DataBim, by.x=c(1,2), by.y=c(1,4))
#  V1       V2  V3.x       V2.y    V3.y V5 V6
#1 15 45420718 black  rs1365242 65.3342  C  A
#write.snp.id.col<-aggregate(V2.y~V3.x, DataCol2,function(x)return(x))
write.snp.id.col<-sapply(unique(DataCol2$V3.x), function(x)return(DataCol2$V2.y[DataCol2$V3.x==x]))
names(write.snp.id.col)<-unique(DataCol2$V3.x)
print(write.snp.id.col)
head(DataCol2)
}else if(!is.null(opt[['pos_ref']])){
write.snp.id.col=list('red'=DataBim$V2[DataBim$V4==opt[['pos_ref']]])
}
pdf(out)
LD.plot2(m, DataBim$V4, draw.chr=TRUE, write.snp.id=TRUE, write.snp.id.col=write.snp.id.col)
dev.off()
