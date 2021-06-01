/*
 *  REMOVE DUPLICATED SNV POSITIONS AND POOR SAMPLES
 *  ================================================
 *
 *  In this step we perfom *pre-QC* filtering of the dataset...
 *
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/intensityPlot.nf"

include {
    getInputChannels;
    checkInputParams;
    getListOfDuplicatePositions;
    removeSamplesWithPoorClinicalData;
    removeDuplicatedSnvPositions;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/duplicatedPositions.nf"

workflow {

    checkInputParams()

    (cohortData,
     samplesWithPoorClinicalData) \
        = getInputChannels()

    duplicatePositions \
        = getListOfDuplicatePositions(
            cohortData)

    sampleFilteredCohortData \
        = removeSamplesWithPoorClinicalData(
            cohortData,
            samplesWithPoorClinicalData)

    removeDuplicatedSnvPositions(
        sampleFilteredCohortData,
        duplicatePositions)
    
}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}