snp.sq <- function(i, j, ld, write.ld.val, color, polygon.par, cex.ld) {
  if(i == j) return(); # pas de plot d'un SNP avec lui même
  if(j < i) { tmp <- j; j <- i; i <- tmp } # i < j
  d <- (j-i)/2
  cx <- i+d
  cy <- -d
  do.call(polygon, c(list(x = cx+c(-1,0,1,0)/2, y = cy + c(0,1,0,-1)/2, col = color), polygon.par))
  if(write.ld.val) text(cx, cy, ld, cex = cex.ld)
}

LD.plot <- function(LD, snp.positions, max.dist = Inf, depth = nrow(LD), graphical.par = list(mar = c(0,0,0,0)),
                    cex.ld, cex.snp, polygon.par = list(border = "white"),
                    color.scheme = function(ld) rgb(1,1-abs(ld),1-abs(ld)),
                    write.snp.id = TRUE, write.ld = function(ld) sprintf("%.2f", ld),
                    draw.chr = TRUE,
                    above.space = 1 + 2*write.snp.id + draw.chr, below.space = 1, 
                    pdf.file, finalize.pdf = TRUE) {

  n <- nrow(LD)
  positions <- if(missing(snp.positions)) rep(0,n) else snp.positions
  # dry run to get graph depth
  graph.depth <- 0
  for(i in seq(1,n-1)) 
    for(j in seq(i+1, min(n,i+depth))) {
      if(positions[j] - positions[i] > max.dist) next;
      graph.depth <- max(graph.depth, (j-i)/2)
    }


  if(!missing(pdf.file)) 
    pdf(pdf.file, width = n/2, height = (graph.depth + above.space + below.space)/2)

  do.call(par, graphical.par) 
  write.ld.val <- ifelse(is.null(write.ld), FALSE, TRUE)
  plot( 0,0, xlim=c(0.5,n+0.5), ylim=c(-graph.depth-below.space,above.space), type="n", 
        asp = 1, xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
  ld <- ""

  if(missing(cex.ld) & write.ld.val) 
    cex.ld = 0.7*par("cex") / max(strwidth(sapply(LD, write.ld)))

  for(i in seq(1,n-1)) 
    for(j in seq(i+1, min(n,i+depth))) {
      if(positions[j] - positions[i] > max.dist) next;
      if(write.ld.val) ld <- write.ld(LD[i,j])
      snp.sq(i, j, ld, write.ld.val, color.scheme(LD[i,j]), polygon.par, cex.ld)
    }
  if(write.snp.id) {
    rs.h <- strheight(rownames(LD))
    if(missing(cex.snp)) cex.snp <- 0.25*par("cex")/max(rs.h)
    text( 1:n, 0, rownames(LD), srt = 90, cex = cex.snp, adj = c(0,0.5))
  }
  if(!missing(snp.positions) & draw.chr) {
    if(write.snp.id) 
      a <- max( strwidth(rownames(LD), cex = cex.snp) )
    else
      a <- 0
    pos <- 1.5 + (n-2)*(snp.positions - snp.positions[1])/diff(range(snp.positions))
    segments( 1:n, a+0.25,  pos, a + 1.5 )
    rect(1.5, a+1.5, n-0.5, a+1.75)
    segments( pos, a+1.5,  pos, a + 1.75 )
  }
  if(!missing(pdf.file) & finalize.pdf) dev.off()
}

