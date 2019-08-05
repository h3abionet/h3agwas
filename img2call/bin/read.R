
library(crlmm)
library(ff)
library(readxl)

sample=read_excel("../aux/dna.xlsx")

print("Read sample sheet")
names = paste(sample$"Array Info.S",sample$"Sentrix ID",sep="_")

names=c("202351420042_R01C01",  "202351420042_R02C01", "202351420042_R03C01",    "202351420042_R04C01",  "202351420042_R05C01", "202351420042_R06C01", "202351420042_R08C01", "202351420049_R01C01", "202351420049_R02C01", "202351420049_R03C01", "202351420049_R04C01","202351420049_R05C01","202351420049_R06C01",  "202351420049_R07C01")

annof =read.table("../aux/A3.annot",sep=",")
print("Read annotation file")




sample$"Sample_ID" = sample$"Institute Sample Label"
sample = sample[ , !(names(sample) %in% drops)]

xx =
     genotype.Illumina(sampleSheet=sample,
     path="/spaces/scott/h3agwas/img2call/aux/allidats",
     arrayNames = names,                  
     arrayInfoColNames=list(barcode="Array Info.S", position="Sentrix ID"),
     highDensity=TRUE, sep="_",
     fileExt=list(green="Grn.idat", red="Red.idat"),
     XY=NULL, anno=annof, genome=37, call.method="krlmm", trueCalls=NULL,
     cdfName='nopackage', copynumber=FALSE, batch=NULL,
     saveDate=FALSE, verbose=TRUE)
