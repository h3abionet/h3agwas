#!/usr/bin/env nextflow
/*
 *  REMOVE LOW QUALITY SAMPLES AND SNVS
 *  ===================================
 *
 *  This script performs a detailed sample and snv quality control for
 *  your GWAS study. 
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
