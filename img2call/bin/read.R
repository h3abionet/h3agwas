
library(crlmm)
library(ff)
library(readxl)

sample=read_excel("../data/testbatch.xlsx")

print("Read sample sheet")
names = paste(sample$"Array Info.S",sample$"Sentrix ID",sep="_")



annof =read.table("../aux/A3.annot",sep=",",header=TRUE)
print("Read annotation file")

path="/spaces/scott/h3agwas/img2call/data/all"
anno=annof
arrayNames=names
arrayInfoColNames=list(barcode="Array Info.S", position="Sentrix ID")
sep="_"
fileExt=list(green="Grn.idat", red="Red.idat")

sample$"Sample_ID" = sample$"Institute Sample Label"


xx =
     genotype.Illumina(sampleSheet=sample,
     path="/spaces/scott/h3agwas/img2call/data/all",
     arrayNames = names,                  
     arrayInfoColNames=list(barcode="Array Info.S", position="Sentrix ID"),
     highDensity=TRUE, sep="_",
     fileExt=list(green="Grn.idat", red="Red.idat"),
     XY=NULL, anno=annof, genome=37, call.method="krlmm", trueCalls=NULL,
     cdfName='nopackage', copynumber=FALSE, batch=NULL,
     saveDate=FALSE, verbose=TRUE)
