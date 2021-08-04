/*
 *  DRAW INTENSITY PLOT FOR A SELECTED SNV
 *  ======================================
 *
 *  This script takes as input an illumina gsgt format dataset, and a 
 *  SNV id chosen by the user, and then draws a plot of the intensity 
 *  (X,Y) values across all samples in the dataset for this chosen
 *  SNV.
 *
 *  The gsgt illumina data format is a set of
 *  genotype reports. A genotpye report is an csv file that contains 
 *  a portion of the genotype results for the cohort. 
 *  There are several genotpye reports for a illumina sequencing 
 *  project, and one file does not correspond to any one snv or any 
 *  one sample.
 *
 *  The script first reads in the genotype reports and, for each
 *  report, extracts the lines relevant to our selected snv and 
 *  adds them to a list. There is then a list for each genotype report
 *  and these lists are merged together. Finally, we draw the
 *  intensity plot from this merged data file. 
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/base.nf"

include {
    getInputChannels;
    checkInputParams;
    filterRecordsForChosenSnv;
    mergeGenotypeReports;
    drawXYintensityPlot;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/intensityPlot.nf"

workflow {

    checkInputParams()

    genotypeReports = getInputChannels()

    genotypeReportsOfChosenSnv \
        = filterRecordsForChosenSnv(
            genotypeReports)

    snvGenotypeReport \
        = mergeGenotypeReports(
            genotypeReportsOfChosenSnv.collect())

    intensityPlot \
        = drawXYintensityPlot(
            snvGenotypeReport)
}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}
