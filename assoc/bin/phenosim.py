#!/usr/bin/env python
##
# a tool for generating a phenotype to a simulated genotype
##

import os
import sys
import getopt
import re
import random
from ph_parser import *
from ph_phenotyper import *
from ph_converter import *


class CommandLine:
	def __init__(self):
		self.file = None
		self.outfile= None
		self.input = 'M'
		self.output = 'E'
		self.quant = 1
		self.num_qtn = 1
		self.diploid = 0
		self.variances = [0.5]
		self.maf_range_causal = [0.05,1]
		self.maf_cutoff = 0.05
		self.num_marker = None
		self.num_genotypes = None
		self.epistasis = 0
		self.epi_effect = 0
		self.epi_hapl = 0
		self.qual_wt_affect = None
		self.qual_der_affect = None
		self.qual_het_affect = None
		self.remove_causal = 0
		self.causal_pos = None
		self.missing_prop= 0.0
		self.dominant=1
		self.pene_tab=None
		self.min_cases=None
		
class CommandlineError(Exception): 
    def __init__(self): 
        pass		


def usage():
	usage='''
	
phenosim - a tool for adding a phenotype to simulated genotypes v0.15


-f			name of input file
-i			type of input file ("G" for GENOME, "M" for ms, msHOT and msms) (default: M)
-o			type of output file ("B" for Blossoc, "E" for EMMA, "P" for PLINK, "Q" for QTDT, "T" for Tassel 3.0; comma separated list of these letters is also supported) (default: E)
-q			logical value if quantitaive (1) or qualitative (0) phenotypes should be simulated (default: 1)
-d			logical value for ploidy, indivduals are either haploid (0) or diploid (1) (default: 0)
--outfile		prefix for the output files (default: name of input file)
--causal_pos		predefined position of causal marker, then the marker closest to this position is selected as causal marker (only possible for a single causal marker) (no default)

--prob_wt		(if q=0) probability of the wild type being affected (no default)
--prob_mut		(if q=0) probability of the (homozygous) mutant being affected (no default)
--prob_het		(if q=0 and d=1) probability of the heterozygous being affected (no default)
--pene_tab		(if q=0 and n=2) filename of penetrance table (no default)
--min_cases		(if q=0) minimum number of cases in the final sample (no default)

-n			(if q=1) number of simulated QTNs (default: 1)
-v			(if q=1) proportion of variance explained by the QTNs (if multiple QTNs are simulated, separate the values by commas, use no spaces) (default: 0.5)
--dominant		(if q=1) additive (0) or dominant (1) model used (default: 0)

-e			(if q=1 and n=2) logical value if epistatic interaction between the two QTNs is wanted (1) or not (0) (default: 0)
--epi_eff		(if q=1 and n=2) proportion of variance explained by the epistatic effect (default: 0)

--epi_hapl		(if n=2) logical value indicating if the epistatic causal markers should lie on a common haploblock (1) or not (0) (default: 0)

--maf_r			MAF range for causal markers (upper and lower bound, separated by a comma, no space) (default: 0.05,1.0)
--maf_c			minimum MAF for sampled markers (default: 0.05)
--n_m			number of sampled markers (default: all simulated markers)
--n_gt			number of sampled genotypes (default: all simulated genotypes)
--miss			proportion of missing alleles (default: 0)
-r			logical value if causal markers should be removed from the sample (1) or not (0) (default: 0)

	
'''
	sys.exit(usage)


def parse_commandline(commandline):
	"""parse commandline and return instance of commandline object"""
	try:
		short_options="f:i:o:q:v:e:r:n:d:"
		long_options="causal_pos=:maf_r=:maf_c=:n_m=:n_gt=:epi_eff=:epi_hapl=:prob_wt=:prob_mut=:prob_het=:miss=:outfile=:dominant=:pene_tab=:min_cases=:"
		long_options=long_options.split(':')
		opts, args = getopt.getopt(commandline, short_options,long_options)
	except getopt.GetoptError:
		sys.stderr.write('Wrong options\n')
		sys.exit(2)

	#print opts,args

	cl = CommandLine()
	for option, argument in opts:
		try:
			option = re.sub('\-+','',option)
			if option in ['f']:
				cl.file = argument
			elif option in ['i']:
				cl.input = argument
				if cl.input not in ['M','G']:
					raise ValueError,'wrong input option'
			elif option in ['o']:
				cl.output = argument.split(',')
				for x in cl.output:
					if x not in ['B','E','P','T','Q', 'N']:
						raise ValueError,'wrong output option'
			elif option in ['outfile']:
				cl.outfile = argument
			elif option in ['q']:
				cl.quant = int(argument)
				if cl.quant not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['n']:
				cl.num_qtn = int(argument)
				if cl.num_qtn<1:
					raise ValueError,'wrong option'
			elif option in ['v']:
				cl.variances = map(float,argument.split(','))
				if sum(cl.variances)>=1.0:
					raise ValueError, 'too much variance'
			elif option in ['causal_pos']:
				cl.causal_pos = float(argument)
				if not (0.0<=cl.causal_pos<=1.0):
					raise ValueError, 'invalid number'
			elif option in ['maf_r']:
				cl.maf_range_causal = map(float,argument.split(','))
				if len(cl.maf_range_causal)!=2:
					raise ValueError, 'invalid range'
			elif option in ['maf_c']:
				cl.maf_cutoff = float(argument)
				if not 0.0<=cl.maf_cutoff<=1.0:
					raise ValueError, 'invalid cutoff'
			elif option in ['n_m']:
				cl.num_marker = int(argument)
				if not (cl.num_marker>0):
					raise ValueError, 'invalid number'
			elif option in ['n_gt']:
				cl.num_genotypes = int(argument)
				if not (cl.num_genotypes>0):
					raise ValueError, 'invalid number'
			elif option in ['e']:
				cl.epistasis = int(argument)
				if cl.epistasis not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['epi_hapl']:
				cl.epi_hapl = int(argument)
				if cl.epi_hapl not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['epi_eff']:
				cl.epi_effect = float(argument)
				if not (0.0<=cl.epi_effect<=1.0):
					raise ValueError, 'invalid number'
			elif option in ['prob_wt']:
				cl.qual_wt_affect = float(argument)
				if not (0.0<=cl.qual_wt_affect<=1.0):
					raise ValueError, 'invalid number'
			elif option in ['prob_mut']:
				cl.qual_der_affect = float(argument)
				if not (0.0<=cl.qual_der_affect<=1.0):
					raise ValueError, 'invalid number'
			elif option in ['prob_het']:
				cl.qual_het_affect = float(argument)
				if not (0.0<=cl.qual_het_affect<=1.0):
					raise ValueError, 'invalid number'
			elif option in ['miss']:
				cl.missing_prop = float(argument)
				if not (0.0<=cl.missing_prop<=1.0):
					raise ValueError, 'invalid number'
			elif option in ['r']:
				cl.remove_causal = int(argument)
				if cl.remove_causal not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['d']:
				cl.diploid = int(argument)
				if cl.diploid not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['dominant']:
				cl.dominant = int(argument)
				if cl.dominant not in [1,0]:
					raise ValueError,'wrong option'
			elif option in ['pene_tab']:
				cl.pene_tab = str(argument)
			elif option in ['min_cases']:
				cl.min_cases = int(argument)
				if cl.min_cases <= 0:
					raise ValueError,'wrong number of cases'
			# include more code here to parse other options
		except ValueError, e:
			err_str=' '.join(['Invalid option:','-'+option,argument])
			raise ValueError, err_str
	return cl

def sample(geno,pos,phen,indices,remove_causal=0,num_snps=0,num_genotypes=0,maf_cutoff=0,missing_prop=0,min_cases=0):

	if min_cases:
		if min_cases>phen.count('1'):
			raise RuntimeError, 'Less cases than minimum number observed, please repeat the coalescence simulation or modify the MAF settings!'

	if remove_causal:

		pos2=[]
		geno2=[]

		for k in xrange(len(geno)):
			geno2.append([])

		for j in xrange(len(pos)):
			if j not in indices:
				pos2.append(pos[j])
				for k in xrange(len(geno)):

					geno2[k].append(geno[k][j])

	else:
		geno2=geno
		pos2=pos


	if num_snps: 
		index_all=range(len(pos2))

		maf_indices=[]

		for i in index_all:
			maf=0.0
			for g in geno:
				if g[i]=='1':
					maf+=1.0/len(geno)
			if maf>=maf_cutoff:
				maf_indices.append(i)
				

	

		if len(maf_indices)<num_snps:
			raise 'not enough SNPs'
						
		#print len(index_all),num_snps
		if remove_causal:
			selected_indices=random.sample(maf_indices,num_snps)
		else:
                        print maf_indices, num_indices
			selected_indices=random.sample(maf_indices,num_snps-len(indices))+indices
		
		selected_indices.sort()
		pos3=[]
		geno3=[]

		for k in xrange(len(geno2)):
			geno3.append([])

		for j in selected_indices:
			#if j not in indices:
			pos3.append(pos2[j])
			for k in xrange(len(geno2)):

				geno3[k].append(geno2[k][j])		
	else:
		pos3=pos2
		geno3=geno2

	if num_genotypes or min_cases: 
	
		if not num_genotypes:
			num_genotypes=len(phen)

		case_samples=[]

		if min_cases:
			case_indices=[]
			for i in xrange(len(phen)):
				if phen[i]=='1':
					case_indices.append(i)
			case_samples=random.sample(case_indices,min_cases)

		index_all=range(len(geno3))

		
		if min_cases:
			for i in case_samples:
				index_all.remove(i)
		else:
			min_cases=0

		
		selected_indices=random.sample(index_all,num_genotypes-min_cases)
		selected_indices=selected_indices+case_samples
		selected_indices.sort()		

		pos4=pos3
		geno4=[]
		phen4=[]

		for k in selected_indices:
			phen4.append(phen[k])

		for k in xrange(num_genotypes):
			geno4.append([])

		for j in xrange(len(pos3)):

			for k in xrange(len(selected_indices)):

				geno4[k].append(geno3[selected_indices[k]][j])		
	else:
		pos4=pos3
		geno4=geno3
		phen4=phen

	num_genotypes=len(geno4)
	num_marker=len(pos4)
	num_alleles=num_marker*num_genotypes
	missing_alleles=num_alleles*missing_prop

	if missing_alleles>0:
		missing_indices=random.sample(range(num_alleles),int(missing_alleles))
		for k in missing_indices:
			i=k/num_genotypes
			j=k%num_genotypes
			geno4[j][i]='NA'




	return  geno4,phen4,pos4

def binary_search(liste,eingabe):
	maxindex = len(liste) - 1
	suche_erfolgreich = False
	index = 0
	while not suche_erfolgreich and index <= maxindex and index >=0:
		mitte = index + ((maxindex - index)/2)
	
		if mitte==maxindex:
			if liste[mitte]>=eingabe:
				return mitte
		if liste[mitte]<=eingabe<liste[mitte+1]:
			return mitte+1
		if mitte==0 and liste[0]>=eingabe:
			return 0
		elif eingabe < liste[mitte]:
			maxindex = mitte - 1
		else:
			index = mitte + 1


def run(cl):

	print 'Thank you for using phenosim!'

	#parse input
	if cl.input=='M':
		genotypes_all,positions_all,raw_all=parse_ms(cl.file,cl.diploid)
	elif cl.input=='G':
		genotypes_all,positions_all,raw_all=parse_genome(cl.file,cl.diploid)
	else:
		raise ValueError

	print 'Number of Simulations: ',len(genotypes_all)

	if 'B' in cl.output and cl.missing_prop>0.0:
		raise ValueError, 'Blossoc supports no missing data'
		
	if ('E' in cl.output or 'B' in cl.output) and cl.diploid:
		raise ValueError, 'only haploids supported for EMMA and BLOSSOC'
		
	if cl.diploid and not cl.quant and not cl.qual_het_affect and cl.num_qtn==1:
		raise ValueError, 'When simulating qualitative phenotypes in diploids, --prob_het has to be defined'
		
	if not cl.quant and cl.num_qtn==2 and not cl.pene_tab:
		raise ValueError, 'No penetrance table given'
		
	if cl.num_qtn!=1 and cl.causal_pos:
		raise ValueError, 'User defined position of QTN only possible for single QTNs'
		
	if not cl.file:
		raise ValueError, 'No input file (-f) defined'


	#now per simulation
	for i in xrange(len(genotypes_all)):

		geno=genotypes_all[i]
                print "geno : ", len(geno)
		pos=positions_all[i]


		#create phenotypes

		causal_index=None

		if cl.causal_pos:
			causal_index=binary_search(pos,cl.causal_pos)

		if cl.quant:

			if not cl.epistasis and (cl.epi_hapl or cl.epi_effect):
				raise ValueError, 'Epistasis has to be selected, when using --epi_hapl or --epi_effect'
				

			if len(cl.variances)!=cl.num_qtn:
				raise ValueError, 'Number of QTNs does not match the number of variances'
                        ### Choix indices

			indices,phen,mafs = assign_phenotype_quant(geno,variance_proportions=cl.variances,epistasis=cl.epistasis,epi_effect=cl.epi_effect,epi_hapl=cl.epi_hapl,maf_range=cl.maf_range_causal,causal_index_pre=causal_index,diploid=cl.diploid,raw=raw_all[i]) 		
		else:

			if (not cl.qual_wt_affect or not cl.qual_der_affect) and cl.num_qtn==1:
				raise ValueError, 'When simulating qualitative phenotypes, probabilities of being affected have to be defined'
				

			if cl.epistasis:
				raise ValueError, 'When simulating qualitative phenotypes, interactions have to be defined in the penetrance table'
				
			if cl.num_qtn==1:

				indices,phen,mafs = assign_phenotype_qual(geno,cl.qual_wt_affect,cl.qual_der_affect,het_affect=cl.qual_het_affect,maf_range=cl.maf_range_causal,causal_index_pre=causal_index,diploid=cl.diploid) 
				
			elif cl.num_qtn==2:

				if not cl.pene_tab:
					raise ValueError, 'No penetrance table given'
			
				f=open(cl.pene_tab)
				lines=f.readlines()
				f.close()
				
				if cl.diploid:
					gts=['0','1','2']
				else:
					gts=['0','1']
					
				if len(lines)!=len(gts):
					raise ValueError, 'Wrong format of penetrance table'
					
				penetrances={}
					
				for j in xrange(len(lines)):
					split=map(float,lines[j].split())
					if len(split)!=len(gts):
						raise ValueError, 'Wrong format of penetrance table'
					for k in xrange(len(gts)):
						penetrances[(gts[j],gts[k])]=split[k]
						penetrances[(gts[k],gts[j])]=split[k]

				indices,phen,mafs = assign_phenotype_qual_2locus(geno,penetrances,maf_range=cl.maf_range_causal,diploid=cl.diploid,raw=raw_all[i],epi_hapl=cl.epi_hapl) 	
			
			
			else:
				raise ValueError, 'When simulating qualitative phenotypes, only one or two causal markers are supported'
			

		#sample marker and individuals



		geno4,phen4,pos4=sample(geno,pos,phen,indices,cl.remove_causal,cl.num_marker,cl.num_genotypes,maf_cutoff=cl.maf_cutoff,missing_prop=cl.missing_prop,min_cases=cl.min_cases)

		#write output
		
		if not cl.outfile:
			cl.outfile=cl.file

		fname='%s%d' %(cl.outfile,i)
		

		if 'E' in cl.output:
			convert2emma(geno4,pos4,phen4,fname)
		if 'B' in cl.output:
			convert2blossoc(geno4,pos4,phen4,fname,diploid=cl.diploid)
		if 'P' in cl.output:
			convert2plink(geno4,pos4,phen4,fname,diploid=cl.diploid)
		if 'T' in cl.output:
			convert2tassel(geno4,pos4,phen4,fname,diploid=cl.diploid)
		if 'Q' in cl.output:
			convert2qtdt(geno4,pos4,phen4,fname,diploid=cl.diploid,quant=cl.quant)
                if 'N' in cl.output :
			PrintJustPheno(geno4,pos4,phen4,fname,diploid=cl.diploid)

		#write causal_file
		f_causal=open('%s.causal' %fname,'w')

		for j in xrange(len(mafs)):
			if cl.quant:
				l='%f\t%d\t%f\t%f\n' %(pos[indices[j]],indices[j],mafs[j],cl.variances[j])

			else:
				l='%f\t%d\t%f\n' %(pos[indices[j]],indices[j],mafs[j])
			f_causal.write(l)

		f_causal.close()




def main():
	if len(sys.argv)==1:
		usage()
	cl = parse_commandline(sys.argv[1:])
	run(cl)


if __name__ == '__main__':
	main()
