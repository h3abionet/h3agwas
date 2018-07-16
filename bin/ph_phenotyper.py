import random
import math
import datetime


def four_gamete_test(marker1,marker2,missing_char='N'):
	if len(marker1)!=len(marker2):
		raise ValueError, 'unequal number of genotypes'

	gametes=[]

	for i in xrange(len(marker1)):
		if marker1[i]==missing_char or marker2[i]==missing_char:
			continue
		if (marker1[i],marker2[i]) not in gametes:
			gametes.append((marker1[i],marker2[i]))


	if len(gametes)<=3:
		return True
	else:
		return False



def assign_phenotype_quant(marker,variance_proportions=[0.7],epistasis=0,epi_effect=0.0,maf_range=[0,1],epi_hapl=1,causal_index_pre=0,diploid=0,dominant=0,raw=[]):




	if epistasis and len(variance_proportions)!=2:
		raise ValueError, 'Epistasis only for two QTNs defined'

	if causal_index_pre and len(variance_proportions)!=1:
		raise ValueError, 'Only one QTN with predefined position possible'

	maf=0.0

	if causal_index_pre:
		causal_index=causal_index_pre

		mafs=[]
		causal_indices=[]


		for geno in marker:
			if geno[causal_index]=='1':
				maf+=1.0/len(marker)

			
		mafs.append(maf)
		causal_indices.append(causal_index)		
			
	
	i=0
	while i<len(variance_proportions) and not causal_index_pre:
		if i==0:			
			causal_indices=[]
			mafs=[]		
		
		i+=1
		maf=0.0

		count=0

		while maf<=maf_range[0] or maf>=maf_range[1]:       
			count+=1

			if count>=2*len(marker[0]):
				raise RuntimeError, 'No SNP in MAF range found'
			
			maf=0.0
			epi_freq=0.0

			
			if epistasis and i>1 and epi_hapl:
				snp1=[]

				
				if diploid:
					marker2=raw
				else:
					marker2=marker
				
				for k in xrange(len(marker2)):
	
					snp1.append(marker2[k][causal_indices[0]])

				end_haploblock=causal_indices[0]+1
					
				for j in xrange(causal_indices[0]+1,int(0.9*len(marker2[0]))):
					snp2=[]
					for k in xrange(len(marker2)):
						snp2.append(marker2[k][j])
					#print j, four_gamete_test(snp1,snp2)
					if not four_gamete_test(snp1,snp2):
						#print j
						end_haploblock=j-1
						break
					else:
						end_haploblock=j
				if causal_indices[0]+2>=end_haploblock:
					i=0
					break

				causal_index=random.randint(causal_indices[0]+1,end_haploblock)

				epi_freq=0.0
				for geno in marker:
					if diploid:
						if geno[causal_indices[0]]!='0' and geno[causal_index]!='0':
					      		epi_freq+=1.0/len(marker)
					else:
						if geno[causal_indices[0]]=='1' and geno[causal_index]=='1':
					      		epi_freq+=1.0/len(marker)
				if epi_freq<=0.000001 or epi_freq>=0.999999:
					
					i=0
					break

			elif epistasis and i>1:
			
				#causal_index=random.randint(int(0.1*len(marker[0])),int(0.9*len(marker[0])))	
				causal_index=random.randint(int((i-1)/float(len(variance_proportions))*len(marker[0])),int((i)/float(len(variance_proportions))*len(marker[0])-1))

				epi_freq=0.0
				for geno in marker:
					if diploid:
						if geno[causal_indices[0]]!='0' and geno[causal_index]!='0':
					      		epi_freq+=1.0/len(marker)
					else:
						if geno[causal_indices[0]]=='1' and geno[causal_index]=='1':
					      		epi_freq+=1.0/len(marker)

  
				if epi_freq<=0.000001 or epi_freq>=0.999999:
					
					i=0
					break
			
						
			else:		
			
				#causal_index=random.randint(int(0.1*len(marker[0])),int(0.9*len(marker[0])))
				causal_index=random.randint(int((i-1)/float(len(variance_proportions))*len(marker[0])),int((i)/float(len(variance_proportions))*len(marker[0])-1))
				print causal_index, len(marker[0]),i

			if causal_index in causal_indices:
				continue

			for geno in marker:
				if diploid:
					if geno[causal_index]=='1':
						maf+=0.5/len(marker)
					elif geno[causal_index]=='2':
						maf+=1.0/len(marker)
				else:
					if geno[causal_index]=='1':
						maf+=1.0/len(marker)

			
		mafs.append(maf)
		causal_indices.append(causal_index)
		#print causal_indices

	phenotypes=[]   
	



	if (math.sqrt(1-sum(map(abs,variance_proportions))-epi_effect))<0:
		raise ValueError,'Illegal Variance Proportions'
	

	for geno in marker:             

		alleles=[]

		for i in causal_indices:
			alleles.append(geno[i])

		alleles=map(int,alleles)
		
		if diploid:
			alleles=map(lambda x:x-1,alleles)

		
			
		phen=math.sqrt(1-sum(map(abs,variance_proportions)))*random.gauss(0,1) 

		for i in xrange(len(alleles)):
			if diploid and not dominant:
				root=variance_proportions[i]/(2*mafs[i]*(1-mafs[i]))
				phen+=alleles[i]*math.sqrt(root)
			else:	
				root=variance_proportions[i]/(mafs[i]*(1-mafs[i]))
				if alleles[i]>0:
					phen+=1*math.sqrt(root)
				else:
					phen+=0*math.sqrt(root)
			

		if epistasis:
		
			if diploid:
				if alleles[0]>0 and alleles[1]>0:	
					root=epi_effect/(epi_freq*(1-epi_freq))
					phen+=1*math.sqrt(root)
			else:		
				if alleles[0] and alleles[1]:
	
					root=epi_effect/(epi_freq*(1-epi_freq))
					phen+=1*math.sqrt(root)
		
		phenotypes.append(phen)

	return causal_indices,phenotypes,mafs    

def assign_phenotype_qual(marker,wt_affect,mut_affect,het_affect=0,maf_range=[0,1],causal_index_pre=0,diploid=0):


	maf=0.0

	if causal_index_pre:
		maf=0.0
		causal_index=causal_index_pre

		for geno in marker:
			if diploid:
				if geno[causal_index]=='1':
					maf+=0.5/len(marker)
				elif geno[causal_index]=='2':
					maf+=1.0/len(marker)						
			elif geno[causal_index]=='1':
				maf+=1.0/len(marker)
		
	count=0
	while (maf<=maf_range[0] or maf>=maf_range[1]) and not causal_index_pre :        

		count+=1
		if count>=len(marker[0]):
			raise RuntimeError, 'No SNP in MAF range found'

		maf=0.0
		causal_index=random.randint(int(0.1*len(marker[0])),int(0.9*len(marker[0])))	

		for geno in marker:
			if diploid:
				if geno[causal_index]=='1':
					maf+=0.5/len(marker)
				elif geno[causal_index]=='2':
					maf+=1.0/len(marker)						
			elif geno[causal_index]=='1':
				maf+=1.0/len(marker)

	phenotypes=[] 

	for geno in marker:   
		if diploid:
			if geno[causal_index]=='0' and random.random()<wt_affect:
				phen='1'			
			elif geno[causal_index]=='1' and random.random()<het_affect:
				phen='1'
			elif geno[causal_index]=='2' and random.random()<mut_affect:
				phen='1'
			else:
				phen='0'
		else:           
			if geno[causal_index]=='0' and random.random()<wt_affect:
				phen='1'			
			elif geno[causal_index]=='1' and random.random()<mut_affect:
				phen='1'
			else:
				phen='0'
		
		phenotypes.append(phen)

	return [causal_index],phenotypes,[maf]   



def assign_phenotype_qual_2locus(marker,penetrances,diploid=0,maf_range=[0,1],epi_hapl=0,raw=[]):


	causal_indices=[]
	mafs=[]
	
	for i in range(2):
	
		maf=0.0
		count=0
	
		while (maf<=maf_range[0] or maf>=maf_range[1]):      
		
			if epi_hapl and i==1:
				snp1=[]

				
				if diploid:
					marker2=raw
				else:
					marker2=marker
				
				for k in xrange(len(marker2)):
	
					snp1.append(marker2[k][causal_indices[0]])

				end_haploblock=causal_indices[0]+1
					
				for j in xrange(causal_indices[0]+1,int(0.9*len(marker2[0]))):
					snp2=[]
					for k in xrange(len(marker2)):
						snp2.append(marker2[k][j])
					#print j, four_gamete_test(snp1,snp2)
					if not four_gamete_test(snp1,snp2):
						#print j
						end_haploblock=j-1
						break
					else:
						end_haploblock=j
				if causal_indices[0]+2>=end_haploblock:
					i=0
					break

				causal_index=random.randint(causal_indices[0]+1,end_haploblock)


			else:

				count+=1

				if count>=len(marker[0]):
					raise RuntimeError, 'No SNP in MAF range found'

				
				causal_index=random.randint(int(0.1*len(marker[0])),int(0.9*len(marker[0])))
				if causal_index in causal_indices:
					continue
					
			maf=0.0
			for geno in marker:
				if diploid:
					if geno[causal_index]=='1':
						maf+=0.5/len(marker)
					elif geno[causal_index]=='2':
						maf+=1.0/len(marker)						
				elif geno[causal_index]=='1':
					maf+=1.0/len(marker)
						
		causal_indices.append(causal_index)
		mafs.append(maf)

	phenotypes=[] 

	for geno in marker:
		pen=penetrances[(geno[causal_indices[0]],geno[causal_indices[1]])]

		
		if random.random()<pen:
			phen='1' 
		else:
			phen='0'
		
		phenotypes.append(phen)


	return causal_indices,phenotypes,mafs   		





