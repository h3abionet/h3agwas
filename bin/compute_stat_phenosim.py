#!/usr/bin/env python3

import sys
import argparse


''' 
compute statistic with p.value output of significance in function of simulate position 
prog --stat file --header_pval headerpvalue --pos_simul file posimul --out file --windows_size[default 10**6 pb]

'''
def read_bim(bim,rs_pos,chro_pos,pos_pos,posa_pos):
    def getpos(x):
        return x[3]
    lire=open(bim)
    liste_info={}
    for ligne in lire :
       spl=ligne.split()
       if spl[chro_pos] not in liste_info:
          liste_info[spl[chro_pos]]=[]
       liste_info[spl[chro_pos]].append([spl[rs_pos],spl[chro_pos],float(spl[pos_pos]),float(spl[posa_pos])])  
    for chro in liste_info :
       liste_info[spl[chro_pos]]= sorted(liste_info[spl[chro_pos]],key=getpos)
    lire.close()
    return liste_info
def get_posinheader(splheader, header, filename, arg):
    if header not in splheader:
       print("col : "+header+" not foud for arg "+arg+" in file "+filename)
       sys.exit(3)
    return splheader.index(header)

## need to be by chro
## return lines

def parseArguments():
    parser = argparse.ArgumentParser(description='compute statistic with p.value output of significance in function of simulate position ')
    parser.add_argument('--stat',type=str, required=True, help="stat file output of gwas")
    parser.add_argument('--header_pval',type=str, required=True, help="pvalue header for stat file")
    parser.add_argument('--header_chro',type=str, required=True, help="chro header for stat file")
    parser.add_argument('--header_pos',type=str, required=True, help="pos header for stat file")
    ## file with no header contains bim information
    parser.add_argument('--pos_simul',type=str, required=True, help="file with list of position simulated, must be in bim format")
    parser.add_argument('--out',type=str, required=True, help="outfile ")
    parser.add_argument('--windows_size',type=str,default="1000000bp", help="windows size aroud position in pos_simul")
    parser.add_argument('--bim',type=str, required=True , help="bim files")
    parser.add_argument('--alpha_lim',type=str, required=True, help="list of alpha used : separated by a comma")
    args = parser.parse_args()
    return args

args=parseArguments()
#4	snp-known111569696	54.21	35954338	A	C	0.245437	9386	0.900235	0.050000
if "bp" in args.windows_size.lower() :
   distance_to_pos=int(args.windows_size.replace("bp",""))
   pos_bim=3
elif "cm" in args.windows_size.lower() :
   distance_to_pos=float(args.around_sig.replace("cm","")) 
   pos_bim=2
else :
   print("args windows_size, unknown : " + args.windows_size+" xbp for baire pair or xcm for centimorgan")
   sys.exit(12)

### read position was simulated
#def read_bim(bim,rs_pos,chro_pos,pos_pos,posa_pos):
pos_sim=read_bim(args.pos_simul,1,0,3,pos_bim)
## append range in pos_sim
for chro in pos_sim :
   for pos in pos_sim[chro] :
      pos+=[pos[3]-distance_to_pos , pos[3]+distance_to_pos]
### read position of bim
pos_all=read_bim(args.bim,1,0,3,pos_bim)
pos_all2={}
for chro in pos_all :
  pos_all2[chro]=[x[2] for x in pos_all[chro]]

## first step we used 
## open file stat
lirestat=open(args.stat)
splheader=lirestat.readline().replace("\n","").split()
poschro=get_posinheader(splheader, args.header_chro, args.stat, "--header_chro")
pospos=get_posinheader(splheader, args.header_pos, args.stat, "--header_pos")
pospval=get_posinheader(splheader, args.header_pval, args.stat, "--header_pval")
listepvalue=[float(x) for x in args.alpha_lim.split(",")]
minpvalue=max(listepvalue)
resbychro_sim=[]
resbychro_nosim=[]
for lines in lirestat :
    spl=lines.split()
    pval=float(spl[pospval])
    if pval<=minpvalue :
       chro=spl[poschro]
       pos=int(spl[pospos])
       pos2=pos_all[chro][pos_all2[chro].index(pos)][3]
       BalIsSim=False
       if chro in pos_sim :
          for possim in pos_sim[chro] :
              if pos2>=possim[4] and pos2<=possim[5] :
                 BalIsSim=True
                 BalSamePos=int(possim[3])==int(pos)
                 InfoPosS=chro+"-"+str(possim[3])
                 resbychro_sim.append([chro,int(pos), pos2,pval, InfoPosS, BalSamePos])
       if BalIsSim==False :
          resbychro_nosim.append([chro,int(pos), pos2,pval])
       #find in bim position

## compute stat
ent=""
listres=""
for pvalue in listepvalue:
    ent+="nsig_simall_"+str(pvalue)+"\tnsig_sim_"+str(pvalue)+"\t"+"nsig_simaround_"+str(pvalue)+"\t"+"nsig_nosim_"+str(pvalue)+"\t"
    nbresbychro_simpval=len(set([x[4] for x in resbychro_sim if x[3]<=pvalue ]))
    nbresbychro_simpval2=len(set([x[4] for x in resbychro_sim if x[3]<=pvalue and x[5]]))
    nbresbychro_simpvalall=len([x[4] for x in resbychro_sim if x[3]<=pvalue])
    nbresbychro_nosimpval=len([x[0] for x in resbychro_nosim if x[3]<=pvalue])
    listres+=str(nbresbychro_simpvalall)+"\t"+str(nbresbychro_simpval)+"\t"+str(nbresbychro_simpval2)+"\t"+str(nbresbychro_nosimpval)+"\t"

ent+="nsnp\tnsnpsim"
listres+=str(sum([len(pos_all[x]) for x in pos_all]))+"\t"+str(sum([len(pos_sim[x]) for x in pos_sim]))

writeout=open(args.out,"w")
writeout.write(ent+"\n")
writeout.write(listres+"\n")
writeout.close()

