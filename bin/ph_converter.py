

def convert2emma(marker,positions,phenotypes,fname):

	marker_lines=[]

	for i in xrange(len(positions)):                        #transpose
		marker_lines.append([])

	for i in xrange(len(marker)):
		for j in xrange(len(marker[i])):
			marker_lines[j].append(marker[i][j])

	for i in xrange(len(marker_lines)):
		marker_lines[i]='\t'.join(map(str,marker_lines[i]))+'\n'     

	f=open('%s.emma_geno' %fname,'w')
	f.writelines(marker_lines)
	f.close()

	f=open('%s.emma_pheno' %fname,'w')
	f.write('\t'.join(map(str,phenotypes))+'\n')
	f.close()
		

	return    




def convert2blossoc(marker,positions,phenotype,fname,diploid=0):

	marker_lines=[]

	for i in xrange(len(marker)):
		marker_lines.append(' '.join(map(str,marker[i]))+'\n')

	for i in xrange(len(marker_lines)):
		marker_lines[i]='%s '%str(phenotype[i]) +marker_lines[i]

	f=open('%s.blossoc_geno' %fname,'w')
	f.writelines(marker_lines)
	f.close()

	for i in xrange(len(positions)):
		positions[i]=float(positions[i])

	f=open('%s.blossoc_pos' %fname,'w')
	f.write(' '.join(map(str,positions))+'\n')
	f.close()

	return

def convert2plink(marker,positions,phenotype,fname,diploid=0):

	f=open('%s.ped'%fname,'w')

	for i in xrange(len(marker)):
		l='%s 1 0 0 0 %s ' %(i,phenotype[i])
		s=' '.join(map(str,marker[i]))+'\n'
		#print s
		if diploid:
			s=s.replace('2','A A')
			s=s.replace('1','A T')
			#print s
			s=s.replace('0','T T')
			s=s.replace('NA','0 0')
		else:
			s=s.replace('1','A A')
			#print s
			s=s.replace('0','T T')
			s=s.replace('NA','0 0')
		#print s
		f.write(l+s)

	f.close()

	f=open('%s.map'%fname,'w')

	for i in xrange(len(positions)):
		l='1 snp%s 0 %s\n' %(i,positions[i])
		f.write(l)

	f.close()

	f=open('%s.pheno'%fname,'w')

	for i in xrange(len(phenotype)):
		l='%s 1 %s\n' %(i,phenotype[i])
		f.write(l)

	f.close()
	

	return

def PrintJustPheno(marker,positions,phenotype,fname,diploid=0):
        f=open('%s.pheno'%fname,'w')

        for i in xrange(len(phenotype)):
                l='%s 1 %s\n' %(i,phenotype[i])
                f.write(l)

        f.close()


        return




def convert2qtdt(marker,positions,phenotype,fname,diploid=0,quant=1):

	f=open('%s.ped'%fname,'w')

	for i in xrange(len(marker)):
		if quant:
			l='1 %s 0 0 2 %s ' %(i,phenotype[i])
		else:
			l='1 %s 0 0 2 %s ' %(i,(int(phenotype[i])+1))
		s=' '.join(map(str,marker[i]))+'\n'
		#print s
		if diploid:
			s=s.replace('2','A A')
			s=s.replace('1','A T')
			#print s
			s=s.replace('0','T T')
			s=s.replace('NA','0 0')
		else:
			s=s.replace('1','A A')
			#print s
			s=s.replace('0','T T')
			s=s.replace('NA','0 0')
		#print s
		f.write(l+s)

	f.close()

	f=open('%s.map'%fname,'w')

	for i in xrange(len(positions)):
		l='1 snp%s 0 %s\n' %(i,positions[i])
		f.write(l)

	f.close()
	
	f=open('%s.dat'%fname,'w')

	if quant:
		f.write('T some_trait\n')
	else:
		f.write('A some_disease\n')

	for i in xrange(len(positions)):
		l='M snp%s\n' %(i)
		f.write(l)

	f.close()


	return

def convert2tassel(marker,positions,phenotype,fname,diploid=0):

	f=open('%s.trait'%fname,'w')
	
	f.write('%s\t1\t1\n'%len(phenotype))

	#f.write('<Trait>\ttrait1\n')
	f.write('trait1\n')
	for i in xrange(len(phenotype)):
		f.write('ind%s\t%s\n' %(i,phenotype[i]))
	f.close()


	f=open('%s.poly'%fname,'w')
	f.write('<Marker>')

	for i in xrange(len(positions)):
		f.write('\t1_%s' %positions[i])
	f.write('\n')

	for i in xrange(len(marker)):
		f.write('ind%s'%i)
		for j in xrange(len(marker[i])):
			allele=marker[i][j]
			if not diploid:
				allele=allele.replace('NA','?')
				allele=allele.replace('0','T')
				allele=allele.replace('1','A')
			else:
				allele=allele.replace('NA','?:?')
				allele=allele.replace('0','T:T')
				allele=allele.replace('1','A:T')			
				allele=allele.replace('2','A:A')
			f.write('\t%s'%allele)
		f.write('\n')

	f.close()

	return




