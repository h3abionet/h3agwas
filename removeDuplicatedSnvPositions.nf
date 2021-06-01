/*
 *  REMOVE DUPLICATED SNV POSITIONS
 *  ===============================
 *
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    getInputChannels;
    checkInputParams;
} from "${projectDir}/modules/manhattanQq.nf"

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/intensityPlot.nf"

include {
    getListOfDuplicatePositions;
    removeDuplicatedSnvPositions;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/duplicatedPositions.nf"

workflow {

    checkInputParams()

    cohortData = getInputChannels()

    duplicatePositions = getListOfDuplicatePositions(cohortData)

    removeDuplicatedSnvPositions(
        cohortData,
        duplicatePositions)

}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}