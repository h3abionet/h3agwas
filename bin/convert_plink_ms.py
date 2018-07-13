#!/usr/bin/env python
import sys
def ReadBim(FileBim):
    LireBim=open(FileBim)
    ListeAltRef=[]
    for Ligne in LireBim :
        Spl=Ligne.split()
        ListeAltRef.append([Spl[4],Spl[5]])
    return  ListeAltRef

def GetCodePos(x, Info):
    if x==Info[0]:
       return "0"
    if x==Info[1]:
       return "1"
    return "0"
 
NomFilePed=sys.argv[1]
NomFileBim=sys.argv[2]
NomFileOut=sys.argv[3]
NomFileOutPed=sys.argv[4]


ListeAltRef=ReadBim(NomFileBim)

Lire=open(NomFilePed)
## Nombre de ligne
CmtInd=0
for Ligne in Lire:
   if CmtInd==0:
      NbSnp=len(Ligne.split())-6
   CmtInd+=1
Lire.close()

NbSnp/=2
print CmtInd, NbSnp, len(ListeAltRef)
Ecrire=open(NomFileOut,'w')
EcrireNewPed=open(NomFileOutPed,'w')
ListeNbSnp=range(0,NbSnp)
if NbSnp!=len(ListeAltRef) :
   sys.exit("snp number between ped : "+str(NbSnp)+" and bim file : "+str(len(ListeAltRef))) 

PosChr=" ".join([str(x/float(NbSnp)) for x in ListeNbSnp])
Ecrire.write("//\n"+"segsites:"+str(NbSnp)+"\n"+"positions: "+PosChr+"\n")
Lire=open(NomFilePed)
for Ligne in Lire :
    GenoCh1=""
    GenoCh2=""
    SplLigne=Ligne.split()#[6::]
    BedLines="\t".join(SplLigne[:6])
    SplLigne=SplLigne[6:]
    for CmtSnp in ListeNbSnp: 
       BedLines+="\t"
       if SplLigne[CmtSnp*2] in  ListeAltRef[CmtSnp] :
          BedLines+=SplLigne[CmtSnp*2]
       else :
          BedLines+=ListeAltRef[CmtSnp][0]
       BedLines+=" "
       if SplLigne[CmtSnp*2+1] in  ListeAltRef[CmtSnp] :
          BedLines+=SplLigne[CmtSnp*2+1]
       else :
          BedLines+=ListeAltRef[CmtSnp][0]
       try :
          GenoCh1+=GetCodePos(SplLigne[CmtSnp*2], ListeAltRef[CmtSnp])
          GenoCh2+=GetCodePos(SplLigne[CmtSnp*2+1], ListeAltRef[CmtSnp])
       except :
          sys.exit("Wrong number of SNPs :"+str(CmtSnp)+" "+str(len(SplLigne))+" ")
    Ecrire.write(GenoCh1+"\n"+GenoCh2+"\n") 
    EcrireNewPed.write(BedLines+"\n")
Lire.close()
Ecrire.close()
EcrireNewPed.close()
    

