/*
 *	DRAW INTENSITY PLOT FOR A SNV
 *	=============================
 *
 *	This script will take as input an illumina gsgt dataset, and a 
 *	SNV id, and then draw a plot of the intensity (X,Y) values across
 *	all samples in the dataset for this fixed SNV. 
 *
 ********************************************************************/

nextflow.enable.dsl=2

params.numberOfInputFiles = 10
params.inputFileDir = '/mnt/lustre3p/groups/CBBI1243/SADaCC/ambroise-data/gwas-cameroon/2018-10-28/'
params.inputFilePrefix = 'WTS_H3Africa_Wonkam_2018.04'

params.snvName = '200610-1'

workflow {

	inputGenotypeDataFiles\
		= channel
			.fromPath(
				params.inputFilePrefix + '_gtReport_File-*')
			.view()

} 