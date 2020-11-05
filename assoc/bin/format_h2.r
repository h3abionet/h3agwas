GetBoltH2<-function(File){
if(length(File)==0)return(data.frame(Gen=NA, GenSd=NA, Env=NA, EnvSd=NA))
Data<-readLines(File)
Env<-grep("h2e ", Data, value=T)
if(length(Env)==0)Env<-c(NA,NA) else Env<-as.numeric(gsub(")","",strsplit(strsplit(Env,split=":")[[1]][2], split=" (", fixed=T)[[1]], fixed=T))
Gen<-grep("h2g (", Data, value=T, fixed=T)
if(length(Gen)==0)Gen<-c(NA,NA) else Gen<-as.numeric(gsub(")","",strsplit(strsplit(Gen,split=":")[[1]][2], split=" (", fixed=T)[[1]], fixed=T))
data.frame(Gen=Gen[1], GenSd=Gen[2], Env=Env[1], EnvSd=Env[2])
}

GetGemmaH2<-function(File){
if(length(File)==0)return(data.frame(Gen=NA, GenSd=NA, Env=NA, EnvSd=NA))
Data<-readLines(File)
Env<-grep("h2e ", Data, value=T)
if(length(Env)==0)Env<-c(NA,NA) else Env<-as.numeric(gsub(")","",strsplit(strsplit(Env,split=":")[[1]][2], split=" (", fixed=T)[[1]], fixed=T))
Gen<-grep("h2g (", Data, value=T, fixed=T)
if(length(Gen)==0)Gen<-c(NA,NA) else Gen<-as.numeric(gsub(")","",strsplit(strsplit(Gen,split=":")[[1]][2], split=" (", fixed=T)[[1]], fixed=T))
data.frame(Gen=Gen[1], GenSd=Gen[2], Env=Env[1], EnvSd=Env[2])
}

GetGemmaH2<-function(File){
if(length(File)==0)return(data.frame(Gen=NA, GenSd=NA, Env=NA, EnvSd=NA))
Env<-c(NA,NA)
Env<-c(NA,NA)
Data<-readLines(File)
Expl<-grep("## pve estimates ", Data,value=T)
if(is.na(Expl))Gen=c(NA, NA)
else{
Expl<-as.numeric(strsplit(Expl, " =")[[1]][2])
ExplSe<-grep("## se(pve)", Data,value=T, fixed=T)
ExplSe<-as.numeric(strsplit(ExplSe, " =")[[1]][2])
}
Gen<-c(Expl,ExplSe)
data.frame(Gen=Gen[1], GenSd=Gen[2], Env=Env[1], EnvSd=Env[2])
}

## probleme ?
GetGCTA<-function(File){
if(length(File)==0)return(data.frame(Gen=NA, GenSd=NA, Env=NA, EnvSd=NA))
Env<-c("NA","NA")
Data<-readLines(File)
if(length(grep('Error:', Data))>0)return(data.frame(Gen=NA, GenSd=NA, Env=NA, EnvSd=NA))
Gen<-grep("Sum of V(G)/Vp", Data, fixed=T, value=T)
if(length(Gen)==0){
Gen<-grep("V(G)/Vp", Data, fixed=T, value=T)
}
if(is.na(Gen)){
Gen<-c('NA','NA')
}else{
Gen<-strsplit(Gen, split="\t")[[1]][c(2,3)]
}
data.frame(Gen=as.numeric(Gen[1]), GenSd=as.numeric(Gen[2]), Env=Env[1], EnvSd=Env[2])
}


GetLDSC<-function(File){
if(length(File)==0)return(data.frame(Gen=NA, GenSd=NA, Env=NA, EnvSd=NA))
#Total Observed scale h2: -0.12 (0.0582)
Env<-c(NA,NA)
Data<-readLines(File)
Gen<-grep("Total Observed scale h2: ", Data, value=T, fixed=T)
if(length(Gen)==0)Gen<-c(NA,NA) else Gen<-as.numeric(gsub(")","",strsplit(strsplit(Gen,split=":")[[1]][2], split=" (", fixed=T)[[1]], fixed=T))
data.frame(Gen=Gen[1], GenSd=Gen[2], Env=Env[1], EnvSd=Env[2])
}

ExtractCorBolt<-function(File){
if(length(File)==0)return(data.frame(Gen=NA, GenSd=NA, Env=NA, EnvSd=NA))
dataFile<-readLines(File)
## extract pheno
namePheno<-gsub(' ','',gsub('\\\\','', gsub('--phenoCol[=]','',grep('--phenoCol', dataFile,value=T))))

h2e<-sapply(strsplit(grep('h2e (', dataFile,value=T,fixed=T),split=':', fixed=T), function(x)return(x[2]))
h2eVal<-as.numeric(sapply(strsplit(h2e, split='(',fixed=T),function(x)x[1]))
h2eSdVal<-as.numeric(sapply(strsplit(h2e, split='(',fixed=T),function(x)gsub(")", "",x[2],fixed=T)))
h2eres<-data.frame(Pheno1=namePheno, Pheno2=namePheno, r2Env=h2eVal,r2SdEnv=h2eSdVal)


h2g<-sapply(strsplit(grep('h2g (', dataFile,value=T, fixed=T),split=':', fixed=T), function(x)return(x[2]))
h2gVal<-as.numeric(sapply(strsplit(h2g, split='(',fixed=T),function(x)x[1]))
h2gSdVal<-as.numeric(sapply(strsplit(h2g, split='(',fixed=T),function(x)gsub(")", "",x[2],fixed=T)))
h2gres<-data.frame(Pheno1=namePheno, Pheno2=namePheno, r2Gen=h2gVal, r2SdGen=h2gSdVal)

h2<-merge(h2eres,h2gres, by=c(1,2))

h2ecor<-sapply(strsplit(grep('resid corr', dataFile,value=T),split=':'), function(x)return(x[2]))
h2ecorval<-sapply(strsplit(grep('resid corr', dataFile,value=T),split=':'), function(x)return(x[1]))
h2ecorVal<-as.numeric(sapply(strsplit(h2ecor, split='(',fixed=T),function(x)x[1]))
h2ecorSdVal<-as.numeric(sapply(strsplit(h2ecor, split='(',fixed=T),function(x)gsub(")", "",x[2],fixed=T)))
NamePos<-matrix(as.integer(unlist(strsplit(gsub(")","",sapply(strsplit(h2ecorval,split="(", fixed=T), function(x)return(x[2]))),split=","))), ncol=2, byrow=T)
r2env<-data.frame(Pheno1=namePheno[NamePos[,1]], Pheno2=namePheno[NamePos[,2]], r2Env=h2ecorVal, r2SdEnv=h2ecorSdVal)

h2gcor<-sapply(strsplit(grep('gen corr', dataFile,value=T),split=':'), function(x)return(x[2]))
h2gcorval<-sapply(strsplit(grep('gen corr', dataFile,value=T),split=':'), function(x)return(x[1]))
h2gcorVal<-as.numeric(sapply(strsplit(h2gcor, split='(',fixed=T),function(x)x[1]))
h2gcorSdVal<-as.numeric(sapply(strsplit(h2gcor, split='(',fixed=T),function(x)gsub(")", "",x[2],fixed=T)))
NamePos<-matrix(as.integer(unlist(strsplit(gsub(")","",sapply(strsplit(h2gcorval,split="(", fixed=T), function(x)return(x[2]))),split=","))), ncol=2, byrow=T)
r2gen<-data.frame(Pheno1=namePheno[NamePos[,1]], Pheno1=namePheno[NamePos[,2]], r2Gen=h2gcorVal, r2SdGen=h2gcorSdVal)

r2<-merge(r2env,r2gen, by=c(1,2))

rbind(h2,r2)
}
ExtractCorGCTA<-function(File){
dataFile<-readLines(File)
rGValue<-grep('^rG', dataFile, value=T)
if(length(rGValue)==0)return(data.frame(r2Env=NA,r2SdEnv=NA, r2Gen=NA, r2SdGen=NA))
rGValue<-matrix(unlist(strsplit(rGValue, split='\t')), ncol=3, byrow=T)
return(data.frame(r2Env=NA,r2SdEnv=NA,r2Gen=mean(as.numeric(rGValue[,2])), r2SdGen=mean(as.numeric(rGValue[,3]))))
}
#
args = commandArgs(trailingOnly=TRUE)
Soft=arg[1]
File=arg[2]
Type=arg[3]
Out=arg[4]
Pheno=arg[5]
if(Soft=="Gemma"){
tmp<-GetGemmaH2(File)
}else if(Soft=="Bolt"){
tmp<-GetBoltH2(File)
}else if(Soft=="GCTA"){
tmp<-GetGCTA(File)
}else if(Soft=="LDSC"){
tmp<-GetLDSC(File)
}

tmp$soft<-Soft
tmp$Type<-Type
tmp$Pheno<-Pheno
write.table(tmp, quote=F, row.names=F, col.names=T,file=Out,sep='\t')
