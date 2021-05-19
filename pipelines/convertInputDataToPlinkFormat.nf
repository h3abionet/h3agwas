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

include {
	illumina2lgen
} from "${params.modulesDir}/converting-data.nf"

workflow {

	inpat = "${params.input_dir}/${params.input_pat}"

	array = Channel.fromPath(params.chipdescription)
	report = Channel.fromPath(inpat).ifEmpty { 
		error "No files match the pattern " + inpat
	}

	message = channel.from('getting', 'data!')

	printMessage(message) | view()

	ped_ch = illumina2lgen(report.combine(array))

}