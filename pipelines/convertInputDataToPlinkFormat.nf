/*
 *	CONVERT INPUT DATA TO PLINK FORMAT
 *	==================================
 *
 *	In this script we take input genotype data in the form provided by
 *	Illumina, namely genome studio genotype (gsgt) and convert it to 
 *	the binary form readable by plink:
 *		+ .bed
 *		+ .bim
 *		+ .fam
 *	such that the data can then be passed through the other steps in the
 *	workflow. 
 *
 **********************************************************************/
nextflow.enable.dsl=2

include {
	printMessage
} from "${params.modulesDir}/hello-world.nf"

workflow {

	message = channel.from('getting', 'data!')

	printMessage(message) | view()

}