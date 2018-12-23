#! /usr/bin/env python

from __future__ import print_function

import argparse
import os
import sys
import struct
from  numpy import empty, uint32,fromfile,uint16

# we avoid the use of backslashes to assist in templatising the code for Nextflow
TAB=unichr(9)
EOL=unichr(10)

FID_nSNPsRead    =   1000   
FID_IlluminaID   =    102   
FID_SD           =    103   
FID_Mean         =    104   
FID_NBeads       =    107   
FID_MidBlock     =    200   
FID_RunInfo      =    300   
FID_RedGreen     =    400   
FID_MostlyNull   =    401   
FID_Barcode      =    402   
FID_ChipType     =    403   

FIDs = [FID_nSNPsRead, FID_IlluminaID, FID_SD, FID_Mean, FID_NBeads, FID_MidBlock, FID_RunInfo, FID_RedGreen, FID_MostlyNull, FID_Barcode, FID_ChipType]

def parseArguments():
    parser=argparse.ArgumentParser()
    parser.add_argument('sample', type=str, metavar='samplesheet'),
    parser.add_argument('idat', type=str, metavar='IDATDIR',help="directory where IDAT files can be found"),
    parser.add_argument('manifest', type=str, metavar='MANIFESTFILE',help="file with Illumina manifest"),
    parser.add_argument('--sample-sample-column',dest="sample_col",type=int,default=0,\
                        help="col in sample file where sample ID can be found (number cols from 0)")
    parser.add_argument("--sentrix-barcode-column",dest="barcode_col",type=int,default=3,\
                        help="col in sample file where Sentrix barccode found (number cols from 0)")
    parser.add_argument("--sentrix-position-column",dest="position_col",type=int,default=4,\
                        help="col in sample file where Sentrix barccode found (number cols from 0)")
    parser.add_argument("--sample-delimiter",dest="sample_delimiter",type=str,default=",",\
                        help="what separates entries in the sample file"),
    parser.add_argument("-a","--allow-missing-idat",dest="allow",action="store_true",default=False,\
                        help="if IDAT files are missing, report only, otherwise crash")
    parser.add_argument("-s","--suppress-warnings",dest="suppress",action="store_true",default=False,\
                        help="suppress warnings -- be careful")
    parser.add_argument("--skip-header",action="store_true",dest="has_header",default=False)
    parser.add_argument("-n","--num-threads",dest="num_threads",type=int,default=1,\
                        help="number of threads for parallelism")
    parser.add_argument("-c","--chrom-pos",dest="chrom_pos",action="store_true",default=False,\
                        help="show snp as chromosome-position (default SNP ID as in manifest")
    parser.add_argument("-o","--out",dest="out",type=str,required=True,\
                        help="name of output file")
    args = parser.parse_args()
    return args


class SampleEntry:
    
    def __init__(self,pid,fs):
        self.pid=pid    # person or sample ID
        self.fs = fs    # red and gree files



class SNP:

    def __init__(self,addr_a,addr_b,strand,name,uid,chrom,pos,the_snps):
        self.addr_a   = addr_a
        self.addr_b   = addr_b
        self.strand   = strand
        self.name     = name
        self.uid      = uid
        self.chrom    = chrom
        self.pos      = pos
        self.alleles  = the_snps
        self.a_pos = self.b_pos = None

    def setPos(self,pos,ab): # position of probes in idat file
        if ab == 0:
            self.a_pos = pos
        else:
            self.b_pos = pos


    def chrom_pos(self):
        return (self.chrom,self.pos)

    def __str__(self):
        return("{}:{}".format(self.chrom,self.pos))


def getiDatHash(idat_dir):
    ''' dict of all idat files and thei locations : different projects organise their 
        idat files differently -- some flat, some in a hierarchical space '''
    tree = os.walk(idat_dir)
    hash = {}
    for (d,subds,fs) in tree:
        for f in fs:
            if f.endswith(".idat"):
                hash[f] = os.path.join(d,f)
    return hash


def parseSampleSheet(args):
    #parse the sample file to extract the IDs of the particpants and their corresponding
    #idat files. Print warning or crash if the files don't exst
    with open(args.sample) as mf:
        idats=[]
        for line in mf:
            recs    = line.split(args.sample_delimiter)
            pid     = recs[args.sample_col]
            barcode = recs[args.barcode_col]
            pos     = recs[args.position_col]
            curr_fs = []
            ok= True
            warning = ""
            for colour in ["Red","Grn"]:
                base_file = "{barcode}_{pos}_{colour}.idat".format(barcode=barcode,pos=pos,colour=colour)
                f =  idat_hash[base_file]
                this_ok = os.access(f,os.R_OK)
                if not this_ok: warning=warning+"Warning: file {} does not exist or readable{}".format(f,EOL)
                ok = ok & this_ok
                curr_fs.append(f)
            if not ok:
                if args.allow: 
                    if not args.suppress:
                       sys.stderr(warning+EOL)
                    continue
                else:
                    sys.exit("Missing idat files: "+EOL+warning)
            idats.append(SampleEntry(pid,curr_fs))
    return idats


def colsOfManifest(fnames):
    ''' return the index(base 0) of the column in the  manifest file for the key fields we need'''
    fields = []
    for name in ["IlmnStrand","Name","SNP","AddressA_ID","AddressB_ID","Chr","MapInfo"]:
        fields.append(fnames.index(name))
    return fields

def getManifest(args):
    # Returns a list of all the SNPs plus an index for each probe saying which SNP it belongs go
    snp_manifest = []
    address_index= {}
    with open(args.manifest) as f:
        line=f.readline()
        while line[:6] != "IlmnID":
            line=f.readline()
        fnames = line.split(",")
        cols=colsOfManifest(fnames)
        oldpos=oldchrom=1
        for line in f:
           fields   = line.split(",")
           if "Controls" in fields[0]: break
           try:
              (strand,name,snps,address_a,address_b,chrom,pos)=map(lambda col: fields[col],cols)
              uid = "{}:{}".format(chrom,pos)
              the_snps = snps[1:-1]
              addr_a   = int(address_a)
              addr_b = int(address_b) if address_b else None
              snp_manifest.append(SNP(addr_a,addr_b,strand,name,uid,chrom,pos,the_snps))
           except IndexError:
              if not args.suppress: sys.stderr.write(line)
        snp_manifest.sort(key=SNP.chrom_pos)
        print("Here",len(snp_manifest))
        for i, snp in enumerate(snp_manifest):
            address_index[snp.addr_a]=(i,0)
            if snp.addr_b>=0: address_index[snp.addr_b]=(i,1)
    print(len(address_index.keys()))
    return (snp_manifest,address_index)



def getNum(f,num_bytes=4):
    #j=i+num_bytes
    if num_bytes==2:
        code='H'
    elif num_bytes==4:
        code='L'
    elif num_bytes==8:
        code='Q'
    data = f.read(num_bytes)
    res, = struct.unpack("<%s"%code,data)
    return res

def getVals(fname):
    with open(fname,"rb") as f:
        # read as string
        #data = f.read()
        magic_number=f.read(4)
        if magic_number != "IDAT":
            sys.exit("Not an IDAT file")
        version = getNum(f)
        if  version != 3:
            sys.exit("IDAT version 3 supported only, found {}".format(version))
        #skip
        getNum(f)
        fcount = getNum(f)
        field_val = {}
        for i in range(fcount):
            fcode = getNum(f,2)
            offset= getNum(f,8)
            field_val[fcode]=offset
            #print(fcode,offset)
        f.seek(field_val[FID_nSNPsRead])
        num_markers =  getNum(f)
        f.seek(field_val[FID_Barcode])
        bcode       =  getNum(f)
        f.seek(field_val[FID_IlluminaID])
        iids =  fromfile(f,dtype=uint32,count=num_markers)
        f.seek(field_val[FID_Mean])
        vals =  fromfile(f,dtype=uint16,count=num_markers)
        return (iids,vals)


def probeIndexInData(idatf,smf,aidx):
    print(idatf)
    probe_addr, intensities = getVals(idatf)
    for i, addr in enumerate(probe_addr):
        try:
            snp_pos,ab = aidx[addr]
            smf[snp_pos].setPos(i,ab)
            #print(i,addr,snp_pos,ab)
        except KeyError:
            if not args.suppress:
                sys.stderr.write("Warning: Probe {} not in manifest{}".format(addr,EOL))
    for snp in smf:
        if not snp.a_pos and not args.suppress:
                sys.stderr.write("Warning: SNP {} not in idat file{}".format(addr,EOL))



def getSNPIntensities(data,s_idx,res,smf):
    ''' For each SNP find probe(s) for that SNP and get values
        data -- numpy array where data to be stored
        sample -- whose file we're dealing with (index)
        res    -- array of red, green values by probe
        smf    -- SNP manififest file '''
    for i, snp in enumerate(smf):
        a_idx = snp.a_pos
        if not a_idx: continue  # warning given earlier
        for colour in [0,1]:
            data[s_idx,colour,i] =  res[colour][a_idx]
        
def showHeading(f,idats):
    f.write("SNP{}Coord{}Alleles".format(TAB,TAB,TAB))
    for entry in idats:
        f.write("{}{}".format(TAB,entry.pid))
    f.write(EOL)

def showSNP(f,data,snp_i,snp,address_pos,AB,num):
   if not address_pos:
       return
   if args.chrom_pos:
      f.write(snp.uid+AB)
   else:
      f.write(snp.name)
   f.write("{}{}{}{}".format(TAB,snp_i,TAB,snp.alleles.replace("/","")))
   for sample_i in range(0,num):
      for colour in [0,1]:
        f.write("{}{}".format(TAB,data[sample_i,colour,snp_i]))
   f.write(EOL)


def showIntensities(args,data,smf,idats,num=-1):
    if num==-1: num=len(data.shape[0])
    old_chrom=f=None
    for snp_i,snp in enumerate(smf):
        if snp.chrom != old_chrom:
            if old_chrom: f.close()
            f=open(args.out+"_"+snp.chrom+".csv","w")
            showHeading(f,idats)
        old_chrom=snp.chrom
        showSNP(f,data,snp_i,snp,snp.addr_a,"",num)
        showSNP(f,data,snp_i,snp,snp.addr_b,"B",num)
        if snp.addr_b:
           showSNP
    f.close()


def batchProcessIDATS(args,data,idats,curr_b,smf):
    batch = range(len(idats))[curr_b]
    for i in batch:
        sample = idats[i].fs
        res = map(lambda fn : getVals(fn)[1], sample)
        print(idats[i].pid)
        getSNPIntensities(data,i,res,smf)



def processIDATS(args,idats,smf,aidx):
     n_samples=len(idats)
     n_snps   = len(smf)
     data     = empty((n_samples,2,n_snps), dtype=uint32)
     probeIndexInData(idats[0].fs[0],smf,aidx)
     for batch in range(args.num_threads):
        curr_b = slice(batch,n_samples,args.num_threads)
        batchProcessIDATS(args,data,idats,curr_b,smf)
     showIntensities(args,data,smf,idats,10)



if __name__ == '__main__':
    args       =  parseArguments()
    print("Reading sample sheet")
    idats      =  parseSampleSheet(args) 
    print("Reading manifest")
    (smf,aidx) =  getManifest(args)
    print("Processing idat files")
    processIDATS(args,idats,smf,aidx)
