#!/usr/bin/env nextflow
/*
 *  REMOVE LOW QUALITY SAMPLES
 *  ==========================
 *
 *  This script performs specific filtering of your cohort data for 
 *  samples specifically. That is, we look for all the samples in the 
 *  cohort data that do not pass typical GWAS quality thresholds and 
 *  then cut them from the data. 
 *  
 *  These quality tests are:
 *      + discordant sex information for a sample
 *      + high snv-missingness rate for a sample
 *      + extreme values of recorded heterozygosisty
 *      + highly related samples
 *  
 *  Each sample has a reported sex in the clinical phenotype fam file
 *  and if this is different from that deduced by plink from the 
 *  genotype data alone, this sample is removed. Samples can have a
 *  large number of missing genotype information which can lead to 
 *  poor or false statistical associations, so these samples are cut. 
 *  
 *  We select the problematic samples in parallel for each of the above
 *  tests, combine these into a single list, and then remove all
 *  samples togther from the cohort data at the end of the script, i.e
 *  *after* drawing all plots. Thus the plots represent the data before
 *  the poor samples are cut.
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    printWorkflowExitMessage;
    collectPlotsTogetherAndZip;
} from "${projectDir}/modules/base.nf"

include {
    getInputChannels;
    checkInputParams;
    getPopulationStratificationReports;
    getDiscordantSampleSexInfoReport;
    getSampleBasedMissingnessReport;
    getIdentityByDescentReport;
    getSampleHeterozygosityReport;
    drawPopulationStratificationPlot;
    drawSampleMissingnessHistogram;
    drawMissingnessHeterozygosityPlot;
    selectSamplesWithDiscordantSexInfo;
    selectSamplesWithHighRelatedness;
    selectSamplesFailingMHTest;
    concatenateSampleLists;
    removeLowQualitySamples;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/sampleQualityControl.nf"


workflow {

    checkInputParams()

    cohortData = getInputChannels()

    cohortEigen \
        = getPopulationStratificationReports(
            cohortData)

    (pcaPlot, eigenvaluePlot) \
        = drawPopulationStratificationPlot(
            cohortEigen)

    cohortSexcheck \
        = getDiscordantSampleSexInfoReport(
            cohortData)

    discordantSexSamples \
        = selectSamplesWithDiscordantSexInfo(
            cohortSexcheck)

    cohortImiss \
        = getSampleBasedMissingnessReport(
            cohortData)

    sampleMissingnessHistogram \
        = drawSampleMissingnessHistogram(
            cohortImiss)

    cohortGenome \
        = getIdentityByDescentReport(
            cohortData)

    highRelatednessSamples \
        = selectSamplesWithHighRelatedness(
            cohortGenome,
            cohortImiss)

    cohortHet \
        = getSampleHeterozygosityReport(
            cohortData)

    missingnessHeterozygosityPlot \
        = drawMissingnessHeterozygosityPlot(
            cohortImiss,
            cohortHet)

    failedMHTestSamples \
        = selectSamplesFailingMHTest(
            cohortImiss,
            cohortHet)

    lowQualitySamples \
        = concatenateSampleLists(channel.empty().mix(
            discordantSexSamples,
            highRelatednessSamples,
            failedMHTestSamples).collect())

    plots = channel
        .empty().mix(
            pcaPlot,
            eigenvaluePlot,
            sampleMissingnessHistogram,
            missingnessHeterozygosityPlot)
        .collect()

    collectPlotsTogetherAndZip(
        "sampleFiltering",
        plots)

    sampleFilteredCohortData \
        = removeLowQualitySamples(
            cohortData,
            lowQualitySamples)

}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}
