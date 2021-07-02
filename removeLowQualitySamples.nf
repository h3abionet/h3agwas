#!/usr/bin/env nextflow
/*
 *  REMOVE LOW QUALITY SAMPLES
 *  ==========================
 *
 *  This script performs a detailed sample and snv quality control for
 *  your GWAS study. 
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/base.nf"

include {
    getInputChannels;
    checkInputParams;
    getDiscordantSampleSexInfoReport;
    getSampleBasedMissingnessReport;
    getHardyWeinbergEquilibriumReport;
    getIdentityByDescentReport;
    getSampleHeterozygosityReport;
    drawSampleMissingnessHistogram;
    drawMissingnessHeterozygosityPlot;
    selectSamplesWithDiscordantSexInfo;
    selectSamplesWithHighRelatedness;
    selectSamplesFailingMHTest;
    concatenateSampleLists;
    removeLowQualitySamples;
    collectPlotsTogetherAndZip;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/sampleQualityControl.nf"


workflow {

    checkInputParams()

    cohortData = getInputChannels()

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

    cohortHwe \
        = getHardyWeinbergEquilibriumReport(
            cohortData)

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

    collectPlotsTogetherAndZip(channel.empty().mix(
        sampleMissingnessHistogram,   
        missingnessHeterozygosityPlot).collect())

    sampleFilteredCohortData \
        = removeLowQualitySamples(
            cohortData,
            lowQualitySamples)

}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}