#!/usr/bin/env nextflow
/*
 *  REMOVE LOW QUALITY SNVS
 *  ===================================
 *
 *  This script removes poor snvs from the cohort data. It is a best
 *  practice in GWAS studies to remove the low quality snvs *after*
 *  removing the snvs, so this is what you must do here. 
 *
 *  The tests that snvs in the cohort data must pass are:
 *      + low missingness rate across the cohort samples
 *      + low p-vaues of differential missingness (i.e. low
 *          associations between missingness and case/control status)
 *      + high minor allele frequencies
 *      + in hardy weinberg equilibrium
 *
 *  As for the sample quality control step before, we select all snvs
 *  first, combine them into one list, and then remove all snvs 
 *  together at the end of the script, thus the plots represent the 
 *  data *before* quality cuts. 
 *
 ********************************************************************/

nextflow.enable.dsl=2

include {
    printWorkflowExitMessage;
    collectPlotsTogetherAndZip;
} from "${projectDir}/modules/base.nf"

include {
    checkInputParams;
    getInputChannels;
    getMissingnessReport;
    drawDifferentialMissingnessPlot;
    selectSnvsWithHighDifferentialMissingness;
    getAlleleFrequencyReport;
    drawAlleleFrequencyPlot;
    selectlowMinorAlleleFrequencySnvs;
    selectSnvsWithHighMissingness;
    getHardyWeinbergEquilibriumReport;
    selectControlOnlyHWETests;
    drawHardyWeinbergEquilibriumPlot;
    selectSnvsOutOfHardyWeinbergEquilibrium;
    concatenateSnvLists;
    removeLowQualitySnvs;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/snvQualityControl.nf"


workflow {

    checkInputParams()

    cohortData = getInputChannels()

    cohortMissing \
        = getMissingnessReport(
            cohortData)

    highMissingnessSnvs \
        = selectSnvsWithHighMissingness(
            cohortMissing)

    differentialMissingnessPlot \
        = drawDifferentialMissingnessPlot(
            cohortMissing)

    highDifferentialMissingnessSnvs \
        = selectSnvsWithHighDifferentialMissingness(
            cohortMissing)

    cohortFrq \
        = getAlleleFrequencyReport(
            cohortData)

    alleleFrequencyPlot \
        = drawAlleleFrequencyPlot(
            cohortFrq)

    lowMinorAlleleFrequencySnvs \
        = selectlowMinorAlleleFrequencySnvs(
            cohortFrq)

    cohortHwe \
        = getHardyWeinbergEquilibriumReport(
            cohortData)

    controlOnlyCohortHwe \
        = selectControlOnlyHWETests(
            cohortHwe)

    controlOnlyHwePlot \
        = drawHardyWeinbergEquilibriumPlot(
            controlOnlyCohortHwe)

    outOfHardyWeinbergEquilibriumSnvs \
        = selectSnvsOutOfHardyWeinbergEquilibrium(
            controlOnlyCohortHwe)

    lowQualitySnvs \
        = concatenateSnvLists(channel.empty().mix(
            highMissingnessSnvs,
            highDifferentialMissingnessSnvs,
            lowMinorAlleleFrequencySnvs,
            outOfHardyWeinbergEquilibriumSnvs).collect())

    snvFilteredCohortData \
        = removeLowQualitySnvs(
            cohortData,
            lowQualitySnvs)

    plots = channel
        .empty().mix(
            differentialMissingnessPlot,
            alleleFrequencyPlot,
            controlOnlyHwePlot)
        .collect()

    collectPlotsTogetherAndZip(
        "snvFiltering",
        plots)

}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}
