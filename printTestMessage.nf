/*
 *  PRINT TEST MESSAGE
 *  ==================
 *
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
	getInputChannels;
	printToScreen;
} from "${projectDir}/modules/testMessage.nf"

workflow {

	message = getInputChannels()

	printToScreen(message) | view()

}
