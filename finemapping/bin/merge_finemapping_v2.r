#!/usr/bin/env Rscript
library("optparse")
library(data.table)
PlotRes<-function(datainwork,DataGenes2 , datagwascatchro,FilePdf,paintfile=NULL, listpval=c('bf_caviar', 'log10bf_fm', 'Posterior_Prob_paint'), listheadprint=c('caviar', 'fm', 'paint'), gwascat=NULL){
 ## defined c
 if(paintfile %in% c('NOFILE',"0","00","01","02", "03","04","05"))paintfile=NULL
 ColI="blue";ColInitial=t_col(ColI,70)
 ColGWASCat=t_col("red4",30);ColFineMap='black';ColCred<-t_col('grey',30)
 datainworkplot<-datainwork
 Chro<-unique(datainworkplot$chromosome)
 datainworkplot$col<-ColInitial
 datainworkplot$baldiff<-F
 datainworkplot$col[datainworkplot$gwascat]<-ColGWASCat
 datainworkplot$baldiff[datainworkplot$gwascat]<-T
 datainworkplot$col[datainwork$is_cred]<-ColCred
 datainworkplot$baldiff[datainwork$is_cred]<-T
 datainworkplot$baldiff[datainworkplot$IsSig]<-T
 datainworkplot$col[datainworkplot$IsSig]<-ColFineMap
 datainworkplot$pch<-16
 datainworkplot$pch[datainworkplot$gwascat]<-15

 datainworkplot$bp_cm<-datainworkplot$position/1000000
 xlim<-range(datainworkplot$bp_cm)

 if(!is.null(paintfile))pdf(FilePdf, width = 7, height = 7*1.5) else pdf(FilePdf,width = 7, height = 7)
 if(!is.null(paintfile))layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,4,4), 9, 2, byrow = TRUE)) else layout(matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3), 8, 2, byrow = TRUE))
 par(mar=c(0,5,2,1))
 plot(datainworkplot$bp_cm[!datainworkplot$baldiff],-log10(datainworkplot$p[!datainworkplot$baldiff]), xlab='bp (cm)', ylab='-log10(Pi)', pch=datainworkplot$pch, cex=1.2, bg=datainworkplot$col[!datainworkplot$baldiff], col=datainworkplot$col[!datainworkplot$baldiff], xaxt='n',ylim=range(-log10(datainworkplot$p))*1.1, frame.plot = FALSE, cex.axis=2,cex.lab=2)
 points(datainworkplot$bp_cm[datainworkplot$baldiff],-log10(datainworkplot$p[datainworkplot$baldiff]), pch=datainworkplot$pch, cex=1.5, bg=datainworkplot$col[datainworkplot$baldiff], col=datainworkplot$col[datainworkplot$baldiff])
 if(any(datainworkplot$gwascat))text(datainworkplot$bp_cm[datainworkplot$gwascat], -log10(datainworkplot$p[datainworkplot$gwascat]), datainworkplot$rsid[datainworkplot$gwascat],  srt=45,col=ColGWASCat, cex=0.8,pos=4)
 text(datainworkplot$bp_cm[datainworkplot$IsSig], -log10(datainworkplot$p[datainworkplot$IsSig]), datainworkplot$rsid[datainworkplot$IsSig],  srt=45,col='darkblue', cex=1.2,pos=4)
 if(any(datainworkplot$is_cred))text(datainworkplot$bp_cm[datainworkplot$is_cred], -log10(datainworkplot$p[datainworkplot$is_cred]), datainworkplot$rsid[datainworkplot$is_cred],  srt=45,col='grey', cex=1.2,pos=4)
 par(mar=c(0,5,0,1))
 plot(xlim, c(0,1), type='n',xaxt='n' ,yaxt='n', ylab='')
 if(nrow(datagwascatchro)>0){
  datagwascatchro$col<-'red'
  datagwascatchro$col[!(datagwascatchro[,headbpgc] %in%  datainworkplot$position)]<-'grey'
  arrows(datagwascatchro$bp_cm,0.8,datagwascatchro$bp_cm,1,code=2,lwd=1.2, length=0.1, col=datagwascatchro$col)
 }
 listgen<-unique(DataGenes2$GENE);nbgenes=length(listgen)
 if(nbgenes>0){
  sapply(1:length(listgen),function(x){
  posy=x/(nbgenes+1)*0.8
  geneinf=DataGenes2[DataGenes2$Type=='gene' & DataGenes2$GENE==listgen[x],]
  lines(geneinf[1,c("begin_cm","end_cm")], c(posy, posy), lty=3)
  CDSInfo<-DataGenes2[DataGenes2$Type=='CDS' & DataGenes2$GENE==listgen[x],]
  if(nrow(CDSInfo)>0)rect(CDSInfo[,"begin_cm"], posy-0.1*posy, CDSInfo[,"end_cm"], posy+0.1*posy, col='black')
  text(geneinf[,"begin_cm"]+(geneinf[,"end_cm"]-geneinf[,"begin_cm"])/2, posy*1.2, geneinf[,'GENE'], col='red', cex=0.75)
 })
 }
 listcol=rep(c('red', 'blue', 'green','black'),10)
 nbplot=length(listpval)
 par(mar=c(0,5,0,1))
 xaxt='n'
 if(is.null(paintfile))xaxt='l'
 plot(xlim, c(0,nbplot+1),  yaxt='n',xlab=paste('chr', Chro,' (Mb)', sep=''), ylab='', type='n', xaxt=xaxt)
 sapply(1:length(listpval), function(x){
  headp<-listpval[x]
  datatmp<-datainworkplot[!is.na(datainworkplot[,headp]) & !is.infinite(datainworkplot[,headp]),]
  datatmp$norm<-(datatmp[,headp]-min(datatmp[,headp], na.rm=T))/(max(datatmp[,headp], na.rm=T)-min(datatmp[,headp], na.rm=T))*50
  datatmp$col2<-apply(cbind(color.gradient(datatmp[,headp], colors=c('white', listcol[x])),datatmp$norm),1,function(x)t_col(x[1],(100-as.numeric(x[2]))))
  points(datatmp$bp_cm ,rep(x-0.5, length(datatmp$bp_cm)),col=datatmp$col2, bg=datatmp$col2,pch=22, cex=3)
  mtext(listheadprint[x],2,at=x-0.5, las=1 , cex=0.5)
 })
 if(!is.null(paintfile)){
  par(mar=c(2,5,0,1))
   plot(xlim, range(datainworkplot$PercAnnot,na.rm=T),ylab='% annotation', type='n')
   #points(datainworkplot$bp_cm, datainworkplot$PercAnnot, type='h', col=color.gradient(-log10(datainworkplot$p),colors=c('white','red')),cex=-log10(datainworkplot$p)/max(-log10(datainworkplot$p))+1)
   balisegc<-datainworkplot$baldiff
   balisecred<-datainwork$is_cred
   balise<-balisecred | balisegc
   points(datainworkplot$bp_cm[!balise], datainworkplot$PercAnnot[!balise], type='h', col='grey',cex=1)
   #points(datainworkplot$bp_cm[balise], datainworkplot$PercAnnot[balise], type='h', col=color.gradient(-log10(datainworkplot$p)+1,colors=c('white','black')), cex=2)
   points(datainworkplot$bp_cm[balisecred], datainworkplot$PercAnnot[balisecred], type='h', col='red', cex=3)
   points(datainworkplot$bp_cm[balisegc], datainworkplot$PercAnnot[balisegc], type='h', col='blue', cex=3)
  }
  layout(matrix(c(1), 1, 2, byrow = TRUE))
  par(mar=c(0,0,0,0))
  lnum=15;postxt=4
  plot(1:lnum, 1:lnum, type='n', xaxt='n',  yaxt='n')
  text(1.1,lnum-1.5, "Legends (from top to bottow):", pos=postxt,adj=0)
  text(1.1,lnum-3.5, "- top box: ", pos=postxt,)
  text(1.1,lnum-4.5, "  - gwas catalog, are square point and rsid red " , pos=postxt)
  text(1.1,lnum-5.5, "  - lead snps defined by gcta are points and rsid in black", pos=postxt)
  text(1.1,lnum-6.5, "  - credible set defined by cond / sss finemapping are in grey", pos=postxt)
  text(1.1,lnum-7.5, "- genes / gwas catalog box : ", pos=postxt,)
  text(1.1,lnum-8.5, "  - arrow : red arrow : present in dataset / grey arrow absent in dataset" , pos=postxt)
  text(1.1,lnum-9.5, "  - genes represented by box (exon) and lines (intron)", pos=postxt)
  text(1.1,lnum-10.5, "- distribution of post probability : ", pos=postxt)
  text(1.1,lnum-11.5, "  - dstribution of post probability for finemapinf caviar, paintor and finemap, color variabtion gave higher of probabilites", pos=postxt)
  text(1.1,lnum-12.5, "  - distribution of annotation : in % for each positions, red correspond to lead position, or credible interval ", pos=postxt)
  if(!is.null(gwascat)){
    print('GC')
    gwascat2<-merge(datainworkplot, gwascat,by.x=c("position"), by.y=c(headbpgc), all=T)
    gwascat2<-gwascat2[!is.na(gwascat2$trait),]
    gwascat2$info=paste(gwascat2$trait,'/PubmedId:',gwascat2$pubMedID)
    #gwascat2$info[gwascat2$IsSig]<-paste(gwascat2$info[gwascat2$IsSig],'/')
    gwascat3<-aggregate(info~position, data=gwascat2, paste, collapse=',') 
    balise<-gwascat3$position %in% datainworkplot$position[datainworkplot$IsSig];gwascat3$info[balise]<-paste('(Sig)',gwascat3$info[balise],sep='')
    balise<-gwascat3$position %in% datainworkplot$position[datainworkplot$is_cred];gwascat3$info[balise]<-paste('(Cred)',gwascat3$info[balise],sep='')
    listgc<-paste(gwascat3$position,gwascat3$info,paste='-')
    lnum<-length(listgc)+1
    layout(matrix(c(1), 1, 2, byrow = TRUE))
    par(mar=c(0,0,0,0))
    plot(0:lnum, 0:lnum, type='n', xaxt='n',  yaxt='n')
    text(1.1,lnum-.5, "Gwas catalog :", pos=postxt)
    if(length(listgc)>0){
    text(rep(1.1, length(listgc)),1:length(listgc), listgc, pos=postxt)
    }
  }
  dev.off()
}

t_col <- function(color, percent = 50, name = NULL) {
 ## Get RGB values for named color
 if(is.na(percent)){
 percent=0
}
rgb.val <- col2rgb(color)

## Make new color using input color as base and alpha set by transparency
t.col2 <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255,alpha = (100 - percent) * 255 / 100,names = name)
## Save the color
invisible(t.col2)
}
color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

getpaintor<-function(datainwork, opt){
 listbayes=list()
 if(!is.null(opt[['paintor']])){
  datapaint=read.table(opt[['paintor']], header=T)
  datapaint$num<-1:nrow(datapaint)
  datapaint<-datapaint[,c(3,2)]
  names(datapaint)<-c('num', 'post_paint')
  listhead<-'post_paint'
  datainwork<-merge(datainwork , datapaint,by='num',all=T)
 }else if(!is.null(opt[['listpaintor']])){
  listfile<-read.table(opt[['listpaintor']], stringsAsFactors=F, sep=";")
  listhead<-c()
 for(CmtF in 1:nrow(listfile)){
  Head<-listfile$V1[CmtF]
  listbayes[[Head]]=as.numeric(readLines(listfile$V3[CmtF]))
  datapaint<-read.table(listfile$V2[CmtF], header=T)
  datapaint$num<-1:nrow(datapaint)
  datapaint<-datapaint[,c(3,2)]
  head<-paste('post_paint_',Head,sep='')
  names(datapaint)<-c('num', head)
  listhead<-c(listhead, head)
  datainwork<-merge(datainwork , datapaint,by='num',all=T)
  datainwork$IsSig[datainwork[,head]>0.5]<-T
 }
 }else{
  print('not found paintor')
  q()
 }
 return(list(data=datainwork, bayes=listbayes, listhead=listhead))
}

getfm<-function(datainwork, opt){
 if(!is.null(opt[['finemap']])){
  datafinemap=read.table(opt[['finemap']], header=T)
  datafinemap<-datafinemap[,c('rsid','chromosome','position','prob','log10bf','mean','sd','mean_incl','sd_incl')]
  names(datafinemap)[-c(1:3)]<-paste(names(datafinemap)[-c(1:3)],'_fm',sep='')
  datainwork=merge(datainwork,datafinemap, by=c('rsid','chromosome','position'), all=T)
  datainwork$IsSig[(!is.infinite(datainwork$log10bf_fm)&datainwork$prob_fm>0.5&datainwork$log10bf_fm>=2)]<-T
  listhead<-'log10bf_fm'
 }else if(!is.null(opt[['listfinemap']])){
  listfile<-read.table(opt[['listfinemap']], stringsAsFactors=F)
  listhead<-c()
 for(CmtF in 1:nrow(listfile)){
  if(length(readLines(listfile$V2[CmtF]))>0){
   Head<-listfile$V1[CmtF]
   datafinemap=read.table(listfile$V2[CmtF], header=T)
   datafinemap<-datafinemap[,c('rsid','chromosome','position','prob','log10bf','mean','sd','mean_incl','sd_incl')]
   names(datafinemap)[-c(1:3)]<-paste(names(datafinemap)[-c(1:3)],'_fm','_',Head,sep='')
   datainwork=merge(datainwork,datafinemap, by=c('rsid','chromosome','position'), all=T)
   datainwork$IsSig[(!is.infinite(datainwork$log10bf_fm)&datainwork$prob_fm>0.5&datainwork$log10bf_fm>=2)]<-T
   listhead<-c(listhead,paste('log10bf_fm','_',Head,sep=''))
   }
  }
  }else{
  print('not found paintor')
  q()
 }
return(list(data=datainwork,  listhead=listhead))
}




 
option_list = list(
  make_option(c("--paintor_fileannot"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-p", "--paintor"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-x", "--listpaintor"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-c", "--cojo"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-d", "--datai"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-b", "--caviarbf"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--finemap"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--listfinemap"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("--listfinemap_cred"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-z", "--paintor_annot"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-g", "--gwascat"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-y", "--headbp_gc"), type="character",  
              help="dataset file name", metavar="character"),
  make_option(c("-w", "--headchr_gc"), type="character", 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt", 
              help="out prefix", metavar="character"),
  make_option(c("-l", "--list_genes"), type="character",
              help="list of genes ", metavar="character")
); 
Test=F
if(Test==F){ 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
}else{
#opt=list('datai'='rs7169818.all','caviarbf'='rs7169818_caviarbf', 'listfinemap'='infooutfinemap.rs7169818', listpaintor='infooutpaintor.rs7169818','gwascat'='~/Data/GWASCat/GWASCat_141019_ckd.tsv', 'cojo'='rs7169818.resgcta', 'list_genes'='gencode.v19.genes', 'paintor_fileannot'='rs7169818.annot','paintor_annot'='listannot.paintor.rs7169818', 'o'='test')
opt_parser = OptionParser(option_list=option_list);
CmdL="--out out --listpaintor infopaintor --cojo 1_55496558_cojo.jma.cojo --datai 1_55496558.all --caviarbf 1_55496558_caviarbf.marginal --list_genes gencode.v19.genes --gwascat gwascat_format_all.csv --headbp_gc chromEnd --headchr_gc chrom --listfinemap infofinemap --paintor_fileannot NOFILE"
commandArgsI=strsplit(CmdL, split=" ")[[1]]
opt=parse_args(opt_parser,args=commandArgsI)

}
headbpgc<-opt[['headbp_gc']];headchrgc<-opt[['headchr_gc']]
if(grep('\\.csv$', opt[['gwascat']])==1)datagwascat<-read.csv(opt[['gwascat']]) else datagwascat<-read.table(opt[['gwascat']],sep='\t', header=T)
 
datagwascat$chropos<-paste(datagwascat[,headchrgc],datagwascat[,headbpgc])
datagwascat$bp_cm<-datagwascat[,headbpgc]/1000000

datai=read.table(opt[['datai']], header=T)
datai$chropos<-paste(datai$chromosome, datai$position)
datai$gwascat<-F
datai$gwascat[datai$chropos %in% datagwascat$chropos]<-T
headz<-'z'
headn<-'n'
if('Z' %in% names(datai))headz<-'Z'
if('N' %in% names(datai))headn<-'N'
datai<-datai[,c("rsid","chromosome","position","allele1","allele2","maf","beta","se", "p",headz,'gwascat',headn)]
datai$logpval<- -log10(datai$p)
datai$num<-1:nrow(datai)
datai$IsSig<-F
datainwork=datai

listheadplot<-c()
listheadprint<-c()

datacav=read.table(opt[['caviarbf']], header=F)
names(datacav)<-c("num", "p_caviarbf")
datainwork=merge(datainwork,datacav, all=T, by='num')
datainwork$IsSig[datainwork$p_caviarbf>0.5]<-T
listheadprint=c(listheadprint,'bf')
listheadplot<-c(listheadplot,'p_caviarbf')

## fine mapping
tmpfm<-getfm(datainwork, opt)

datainwork<-tmpfm[['data']]
listheadprint=c(listheadprint , 'fm (sss)')
listheadplot<-c(listheadplot,tmpfm[['listhead']])
##
infopaint<-getpaintor(datainwork, opt)
datainwork<-infopaint[['data']]
listheadprint=c(listheadprint ,'paintor')
listheadplot<-c(listheadplot, infopaint[['listhead']])

## data for  gcta
headgcta=opt[['cojo']]
datagcta2<-read.table(headgcta,header=T)
datagcta2$logpJ<- -log10(datagcta2$pJ)
datainwork<-merge(datainwork,datagcta2[,c('Chr','SNP','bp', 'bJ','bJ_se','pJ','LD_r', 'logpJ')], by=c('rsid','chromosome','position'), by.y=c('SNP','Chr','bp'), all=T)
datainwork$IsSig[!is.na(datainwork$pJ)]<-T


if(!is.null(opt[['paintor_fileannot']]) & !(opt[['paintor_fileannot']] %in% c('NOFILE', "0","00","01"))){
	DataAnnot<-read.table(opt[['paintor_fileannot']], header=T)
        print(!is.null(opt[['paintor_annot']]))
	if(is.null(opt[['paintor_annot']])) HeadAnnot = names(DataAnnot)
        else {
          if(opt[['paintor_annot']] %in% c('NOFILE', "0","00","01","02"))  HeadAnnot = names(DataAnnot)
          else  HeadAnnot<-readLines(opt[['paintor_annot']])
        }
	HeadAnnot2<-gsub("-",".",HeadAnnot)
	DataAnnot$num<-1:nrow(DataAnnot)
	DataAnnot<-merge(datainwork[,c('num','rsid','chromosome','position')], DataAnnot[,c('num',HeadAnnot2)], by='num')
	DataAnnotPerc<-data.frame(num=DataAnnot[,c('num')], PercAnnot=apply(DataAnnot[,-c(1:4)],1, function(x)sum(x)/length(x)*100))
	datainwork=merge(datainwork,DataAnnotPerc, by='num', all=T)
}
datainwork$is_cred<-F
if(!is.null(opt[['listfinemap_cred']])){
 listcred<-read.table(opt[['listfinemap_cred']])
 listcred_head<-c()
 for(cmt in 1:nrow(listcred)){
  infofile<-read.table(as.character(listcred[cmt,2]), header=T,stringsAsFactors=F)
  head<-paste('cred_',listcred[cmt,2],sep='')
  datainwork[,head]<-F
  print(as.vector(as.matrix(infofile)))
  datainwork[datainwork$rsid %in% as.vector(as.matrix(infofile)),head]<-T
  print(table(datainwork[,head]))
  listcred_head<-c(listcred_head,head)
 }
 datainwork$is_cred<-F 
 if(length(listcred_head)>1)datainwork$is_cred<-apply(datainwork[,listcred_head], 1, any)
 else if(length(listcred_head)>0)datainwork$is_cred<-datainwork[,listcred_head]
}


Chro<-as.character(unique(datainwork$chromosome))
xlimcm<-range(datainwork$position)/1000000
## 
DataGenes<-fread(opt[['list_genes']])
DataGenes$begin_cm<-DataGenes$BEGIN/1000000
DataGenes$end_cm<-DataGenes$END/1000000
DataGenes2<-as.data.frame(DataGenes[as.character(DataGenes$CHR)==as.character(Chro) & ((DataGenes$begin_cm>=xlimcm[1] & DataGenes$begin_cm<=xlimcm[2]) | (DataGenes$end_cm>=xlimcm[1] & DataGenes$end_cm<=xlimcm[2]) | (xlimcm[1]>=DataGenes$begin_cm & xlimcm[1]<=DataGenes$end_cm) | (xlimcm[2]>=DataGenes$begin_cm & xlimcm[2]<=DataGenes$end_cm)),])
datagwascatchro<-datagwascat[as.character(datagwascat[,headchrgc])==as.character(Chro) & datagwascat$bp_cm>=xlimcm[1] & datagwascat$bp_cm<=xlimcm[2],]

FilePdf=paste(opt[['out']],'.pdf',sep='')
PlotRes(datainwork,DataGenes2 , datagwascatchro,FilePdf,opt[['paintor_fileannot']], listpval=listheadplot, listheadprint=listheadprint,  gwascat=datagwascatchro)



datainwork$Gene<-sapply(datainwork$position, function(x, DataGenes){
paste(unique(DataGenes[x>=DataGenes$BEGIN & x<=DataGenes$END ,'GENE']), collapse=',')
}, DataGenes2)

datainwork$CDS<-sapply(datainwork$position, function(x, DataGenes){
paste(unique(DataGenes[x>=DataGenes$BEGIN & x<=DataGenes$END ,'GENE']), collapse=',')
}, DataGenes2[DataGenes2$Type=='CDS' , ])
datainwork$CDS[datainwork$CDS=='']<-NA
datainwork$Gene[datainwork$Gene=='']<-NA

if(nrow(datagwascatchro)>0){
datagwascatchro$found<-T
datagwascatchro$found[!(datagwascatchro[,headbpgc] %in%  datainwork$position)]<-F
}
write.table(datainwork, file=paste(opt[['out']],'.all.out', sep='') ,sep='\t', row.names=F, col.names=T, quote=F)
write.table(datagwascatchro, file=paste(opt[['out']],'.gwascat.out', sep='') ,sep='\t', row.names=F, col.names=T, quote=F)
