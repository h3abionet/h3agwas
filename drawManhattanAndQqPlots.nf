/*
 *	DRAW MANHATTAN and QQ PLOTs
 *	===========================
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
	getInputChannels;
	getAssociationReport;
	drawManhattanPlot;
	drawQqPlot;
	sendWorkflowExitEmail;
} from "${projectDir}/modules/manhattanQq.nf"

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/intensityPlot.nf"

workflow {

	cohortData = getInputChannels()

	associationReport = getAssociationReport(cohortData.collect())

	drawManhattanPlot(associationReport)

	drawQqPlot(associationReport)

}

workflow.onComplete {
    printWorkflowExitMessage()
    //sendWorkflowExitEmail()
}