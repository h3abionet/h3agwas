#!/usr/bin/env python3
import sys
import re
args = sys.argv
print(args)
if len(sys.argv)==0 :
  filegenes="/spaces/jeantristan/GWAS/Ressource/gencode.v19.annotation.gtf"
else :
 filegenes=args[1]

readgene=open(sys.argv[1])
writegene=open("gencode.v19.genes", 'w')
writegene.write("\t".join(["CHR", "Type","BEGIN", "END","GENE"])+'\n')
patrn='gene_name'
for line in readgene :
 if line[0] != "#" :
   splline=line.split('\t')  
   #['chr1', 'ENSEMBL', 'CDS', '6186632', '6186806', '.', '-', '0', 'gene_id "ENSG00000116254.13"; transcript_id "ENST00000378021.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "CHD5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "CHD5-201"; exon_number 26; exon_id "ENSE00003551672.1"; level 3; protein_id "ENSP00000367260.1"; tag "basic"; havana_gene "OTTHUMG00000000952.6";\n']^C
   if splline[2]=='CDS' or splline[2]=='gene' :
      tmp=splline[8].split(';')
      resgene=[x for x in tmp if patrn in x ]
      if len(resgene)> 0:
        resgene=resgene[0].replace(patrn, '').replace(' ','').replace('"',"")
        chro=resgene[0].replace('chr','')
        writegene.write("\t".join([chro, splline[2], splline[3], splline[4], resgene])+'\n')

writegene.close()
readgene.close()

