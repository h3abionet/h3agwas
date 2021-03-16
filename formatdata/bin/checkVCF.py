#!/usr/bin/env python
import sys, os
import logging

VERSION = "version 1.4 (20140115)"

# all is a keyword since Python 2.7
try:
    all
except:
    def all(iterable):
        for element in iterable:
            if not element:
                return False
            return True

# convenient functions
def myopen(fn):
    import gzip
    f = gzip.open(fn)
    try:
        f.read(2)
        f.close()
        return gzip.open(fn)
    except:
        f.close()
        return open(fn)
    
def checkGTformat(gt):
    if gt == '.': return True
    genos = gt.replace('|','/').split('/')
    if len(genos) != 1 and len(genos) != 2: return False
    for g in genos:
        if g == '.': continue
        if g.isdigit(): continue
        return False
    return True

#assuming gt is valid, return a integer
def getGeno(gt):
    if gt.find('.') >= 0:
        return -1
    genos = gt.replace('|','/').split('/')
    return sum ( [int(g) for g in genos])

# .fai format
# contig, size, location, basesPerLine, bytesPerLine
# e.g.
# 1       249250621       52      60      61
# 2       243199373       253404903       60      61
class GenomeSequence:
    def __init__(self):
        self.fn = ''
        self.data = {}
        self.fileHandle = None
        # for earlier Python version, we need to define os.SEEK_SET manually
        try:
            os.SEEK_SET
        except AttributeError:
            os.SEEK_SET, os.SEEK_CUR, os.SEEK_END = range(3)
            
    def open(self, fn):
        # read index
        if not os.path.exists(fn + '.fai'):
            print >> sys.stderr, "Cannot .fai index file for %s, consider create index using 'samtools faidx %s'" % (fn, fn)
            return False
        self.fn = fn
        try:
            self.fileHandle = open(self.fn)
            for ln in myopen(fn + '.fai'):
                fd = ln.strip().split()
                if fd[0][:3].upper() == "CHR":
                    fd[0] = fd[0][3:]
                key = fd[0]
                val = map(int, fd[1:])
                self.data[key] = val
                self.data['CHR' + key] = val
            return True
        except:
            return False
    # return base at 0-based position
    # return None if something wrong
    def getBase(self, chrom, pos):
        chrom = chrom.upper()
        if chrom not in self.data:
            return None
        size, loc, basesPerLine, bytePerLine = self.data[chrom]
        if pos < 0 or pos >= size:
            return None
        lineNo = pos / basesPerLine
        remainder = pos % basesPerLine
        filePos = loc + lineNo * bytePerLine + remainder
        #print filePos, lineNo, remainder, self.data[chrom]
        self.fileHandle.seek(filePos, os.SEEK_SET)
        return self.fileHandle.read(1)
    def close(self):
        if self.fileHandle:
            self.fileHandle.close()

def actionItem(logger = sys.stderr):
    print >> logger, "---------------     ACTION ITEM     ---------------"    

def usage():
    print("Usage: ")
    print("%s -r ref.fa -o preifx input.vcf: check VCF for strand" % sys.argv[0] )

class Logger:
    def __init__(self, fn):
        logging.basicConfig(level=logging.DEBUG,
                            #               format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                            format='%(asctime)s %(message)s',               
                            datefmt='%m-%d %H:%M',
                            filename=fn,
                            filemode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        formatter = logging.Formatter('%(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)
    def write(self, msg):
        if msg == '\n': return
        logging.info(msg)

def banner(logger = sys.stderr):
    print >> logger, "checkVCF.py -- check validity of VCF file for meta-analysis"
    print >> logger, VERSION
    print >> logger, "contact zhanxw@umich.edu or dajiang@umich.edu for problems."
    
if __name__ == '__main__':
    try:
        import getopt
        optlist, args = getopt.getopt(sys.argv[1:], 'r:o:')
        optlist = dict(optlist)
        refFile = optlist['-r']
        outPrefix = optlist['-o']
        if len(args) == 1:
            vcfFile = args[0]
        else:
            print >> sys.stderr, "Please provide one VCF at a time"
            sys.exit(1)
    except:
        usage()
        print
        raise
        sys.exit(1)

    # check ref
    # check genotypes
    # check af sites where AF > 0.5

    # from XiaoweiLib import GenomeSequence
    gs = GenomeSequence()
    if not gs.open(refFile):
        print >> sys.stderr, "Cannot open reference genome, exiting..."
        sys.exit(1)

    logger = Logger(outPrefix + '.check.log')
    banner(logger)
    
    fDup = open(outPrefix + '.check.dup', 'wt')
    fRef = open(outPrefix + '.check.ref', 'wt')
    fNonSnp = open(outPrefix + '.check.nonSnp', 'wt')
    fMono = open(outPrefix + '.check.mono', 'wt')
    fGeno = open(outPrefix + '.check.geno', 'wt')
    fAF = open(outPrefix + '.check.af', 'wt')

    nRef, nGeno, nAF = 0, 0, 0
    nNonSnp = 0
    snpSite = set()
    nDupSite = 0
    nMono = 0
    
    nField = -1
    ACGT = set(['A', 'C', 'G', 'T'])
    ACGTM = set(['A', 'C', 'G', 'T', '.'])
    lineNo = -1

    vcfHeaderLine = 0
    vcfSiteLine = 0
    vcfSample = 0
    chrWarningGiven = False
    print >> logger, "Python version is [ %s ] " % '.'.join(map(str, sys.version_info)).strip()
    print >> logger, "Begin checking vcfFile [ %s ]" % vcfFile
    try:
        prevChrom, prevPos = None, None
        for lineNo, ln in enumerate(myopen(vcfFile)):
            if lineNo % 10000 == 0 and lineNo != 0:
                print >> logger, "[ %d ] lines processed \r" % lineNo,
            if not ln or ln.startswith('##'):
                vcfHeaderLine += 1
                continue
            if ln.startswith('#CHROM'):
                vcfHeaderLine += 1
                fd = ln.strip().split()
                nField = len(fd)
                # check duplicated sample ids
                if len(set(fd[9:])) != len(fd[9:]):
                    actionItem(logger)
                    print >> logger, "Your VCF file have duplicated sample IDs, please fix them and re-run checkVCF.py"
                    sys.exit(1)
                vcfSample = len(fd) - 9
                # XX
                ##print '%d sample loaded' % (vcfSample)
                ##print 'nField = %d' % nField
                continue

            if ln[:3].upper() == "CHR" and not chrWarningGiven:
                chrWarningGiven = True
                # sys.exit(1)

            vcfSiteLine += 1
            fd = ln.strip().split()
            if len(fd) != nField:
                print >> logger, "Line [ %d ] does not have correct column number, exiting!" % (lineNo + 1)
                print >> logger, "Current line has %d columns." % (len(fd))
                print >> logger, "First 50 characters in the current line content [ %s ]. " % (ln.strip()[:50])
                sys.exit(1)

            chrom, pos, rsId, ref, alt, qual, filt, info, format = fd[:9]
            if len(ref) != 1 or len(alt) != 1:
                fNonSnp.write("%s\n" % ('\t'.join(fd[:5])))
                nNonSnp += 1
                continue
            if ref not in ACGT or alt not in ACGTM:
                fNonSnp.write("%s\n" % ('\t'.join(fd[:5])))
                nNonSnp += 1        
                continue

            site = '%s:%s' % (chrom, pos)
            if site in snpSite:
                print >> logger, "Duplicated site [ %s ]" % site
                fDup.write('DuplicatedSite\t%s:%s\n' % (chrom, pos))
                nDupSite += 1
                continue
            else:
                snpSite.add(site)

            # check VCF in assending order
            if prevChrom != chrom:
                prevChrom, prevPos = chrom, int(pos)
            else:
                if prevPos > chrom:
                    actionItem(logger)
                    print >> logger, "At line [ %d ], genomic position %s:%s is before previous position %s:%s " (lineNo + 1, chrom, pos, prevChrom, prevPos)
                    continue
            

            # check ref
            trueRef = gs.getBase(chrom, int(pos) - 1)
            if trueRef == None:
                fRef.write('FailedGetBase\t%s:%s\n' % (chrom, pos))
                nRef += 1
                continue
            if ref != trueRef:
                fRef.write('MismatchRefBase\t%s:%s:%s-%s/%s\n' % (chrom, pos, trueRef, ref, alt))
                nRef += 1
                continue

            # check genotype
            try:
                gtIndex = [idx for idx, i in enumerate(format.split(':')) if i == 'GT'][0]
            except:
                print >> logger, "Line [ %d ] does not have GT defined in the FORMAT field"
                continue

            # check genotype
            try:
                genos = [i.split(':')[gtIndex] for i in fd[9:]]
            except:
                fGeno.write('IndividualMissingGTField\tLine:%d\n' % (lineNo + 1) )
                nGeno += 1
                continue
            if not all([ checkGTformat(g) for g in genos]):
                fGeno.write('IndividualHasInvalidGT\tLine:%d\n' % (lineNo + 1) )
                nGeno += 1
                continue
            
            # check AF
            geno = [getGeno(g) for g in genos]
            ac = sum( (g for g in geno if g > 0) )
            nSample = sum( (1 for g in geno if g >= 0 ) )
            if ac == 0 or ac == 2*nSample:
                # monomorphic site
                fMono.write('%s:%s\t%d\t%d\n' % (chrom, pos, ac, nSample))
                nMono += 1

            if nSample > 0:
                af = 1.0 * ac / nSample / 2
            else:
                af = 0.0
            if af > 0.5:
                fAF.write('%s:%s\t%s\t%s\t%f\n' % (chrom, pos, ref, alt, af))
                nAF += 1
    except SystemExit:
        sys.exit(1)
    except KeyboardInterrupt:
        print >> logger, "VCF checking has been stopped at line [ %d ]" % (lineNo + 1)
        print >> logger, " [ %s ... ] " % ln[:50]
        sys.exit(1)
    except IOError:
        print >> logger, "IOError happened..."
        raise
        sys.exit(1)
    except Exception as e:
        print >> logger, "VCF checking failed at line [ %d ]" % (lineNo + 1)
        print >> logger, " [ %s ... ] " % ln[:50]
        print >> logger, "Python exceptions occurred [ %s ]!" % e
        print >> logger, "Please report the above to zhanxw@gmail.com"
        raise
        sys.exit(1)

    if chrWarningGiven:
        print >> logger, "---------------     WARNING     ---------------"
        print >> logger, "Detected that chromosome names have 'chr' prefix..."
        print >> logger, "Please consider using the following command to clean your VCF file and then re-run checkVCF.py"
        print >> logger, '(grep ^"#" $your_old_vcf; grep -v ^"#" $your_old_vcf | sed \'s:^chr::ig\' | sort -k1,1n -k2,2n) | bgzip -c > $your_vcf_file '
        
    print >> logger, "---------------     REPORT     ---------------"
    print >> logger, "Total [ %d ] lines processed" % (lineNo + 1)
    print >> logger, "Examine [ %d ] VCF header lines, [ %d ] variant sites, [ %d ] samples" % (vcfHeaderLine, vcfSiteLine, vcfSample)
    print >> logger, "[ %d ] duplicated sites" % (nDupSite)
    print >> logger, "[ %d ] NonSNP site are outputted to [ %s ]" % (nNonSnp, outPrefix + '.check.nonSnp')
    print >> logger, "[ %d ] Inconsistent reference sites are outputted to [ %s ]" % (nRef, outPrefix + '.check.ref')
    print >> logger, "[ %d ] Variant sites with invalid genotypes are outputted to [ %s ]" % (nGeno, outPrefix + '.check.geno')
    print >> logger, "[ %d ] Alternative allele frequency > 0.5 sites are outputted to [ %s ]" % (nAF, outPrefix + '.check.af')
    print >> logger, "[ %d ] Monomorphic sites are outputted to [ %s ]" % (nMono, outPrefix + '.check.mono')

    fDup.close()
    fRef.close()
    fNonSnp.close()
    fMono.close()
    fGeno.close()
    fAF.close()

    actionItem(logger)
    if nDupSite > 0 or nRef > 0 or nGeno > 0:
        if nDupSite > 0:
            print >> logger, "* Remove duplicated sites and rerun checkVCF.py"
        if nRef > 0:
            print >> logger, "* Read %s.check.ref, for autosomal sites, make sure the you are using the forward strand" % outPrefix
        if nGeno > 0:
            print >> logger, "* Open %s.check.geno, using line number there to examine the original VCF files; make sure genotypes are correct." % outPrefix
    else:
        print >> logger, "* No error found by checkVCF.py, thank you for cleanning VCF file."
    print >> logger, "* Upload these files to the ftp server (so we can double check): %s.check.log %s.check.dup %s.check.noSnp %s.check.ref %s.check.geno %s.check.af %s.check.mono" % ( (outPrefix,)*7)
    
