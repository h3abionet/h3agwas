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

    getInputChannels;
    checkInputParams;
    getListOfDuplicatePositions;
    removeSamplesWithPoorClinicalData;
    removeDuplicatedSnvPositions;
    identifyIndivDiscSexinfo;
    generateSnpMissingnessPlot;
    generateIndivMissingnessPlot;
    getInitMAF;
    showInitMAF;
    showHWEStats;
    removeQCPhase1;
    pruneForIBD;
    findRelatedIndiv;
    calculateSampleHeterozygosity;
    generateMissHetPlot;
    getBadIndivsMissingHet;
    noSampleSheet;
    removeQCIndivs;
    calculateSnpSkewStatus;
    generateDifferentialMissingnessPlot;
    findSnpExtremeDifferentialMissingness;
    removeSkewSnps;
    calculateMaf;
    generateMafPlot;
    findHWEofSNPs;
    generateHwePlot;
    collectPlotsTogetherAndZip;
    sendWorkflowExitEmail;

} from './modules/qualityControl.nf';

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/intensityPlot.nf"


workflow {

    checkInputParams()

    (cohortData, poorSamples, additionalPhenotypes) = getInputChannels()

    DuplicatePositions = getListOfDuplicatePositions(cohortData)

    RemovePoorSamples = removeSamplesWithPoorClinicalData(cohortData,poorSamples)

    RemoveDuplicateSnps = removeDuplicatedSnvPositions(RemovePoorSamples,DuplicatePositions)

    IdentifyIndivDiscSexinfo = identifyIndivDiscSexinfo(RemoveDuplicateSnps)

    GenerateSnpMissingnessPlot = generateSnpMissingnessPlot(IdentifyIndivDiscSexinfo)

    GenerateIndivMissingnessPlot = generateIndivMissingnessPlot(IdentifyIndivDiscSexinfo)

    GetInitMAF = getInitMAF(cohortData)

    ShowInitMAF = showInitMAF(GetInitMAF)

    ShowHWEStats = showHWEStats(IdentifyIndivDiscSexinfo)

    RemoveQCPhase1 = removeQCPhase1(cohortData)

    PruneForIBD  = pruneForIBD(RemoveQCPhase1)

    FindRelatedIndiv = findRelatedIndiv(IdentifyIndivDiscSexinfo,PruneForIBD)

    SampleHeterozygosity = calculateSampleHeterozygosity(RemoveQCPhase1)

    GenerateMissHetPlot = generateMissHetPlot(SampleHeterozygosity)

    BadIndivsMissingHet = getBadIndivsMissingHet(SampleHeterozygosity)

    NoSampleSheet = noSampleSheet()

    RemoveQCIndivs = removeQCIndivs(BadIndivsMissingHet,FindRelatedIndiv,IdentifyIndivDiscSexinfo,
                                    NoSampleSheet, RemoveQCPhase1)

    CalculateSnpSkewStatus = calculateSnpSkewStatus(RemoveQCIndivs,additionalPhenotypes)

    DifferentialMissingnessPlot = generateDifferentialMissingnessPlot(CalculateSnpSkewStatus)

    FindSnpExtremeDifferentialMissingness = findSnpExtremeDifferentialMissingness(CalculateSnpSkewStatus)

    RemoveSkewSnps = removeSkewSnps(RemoveQCIndivs,FindSnpExtremeDifferentialMissingness)

    CalculateMaf = calculateMaf(RemoveSkewSnps)

    GenerateMafPlot = generateMafPlot(CalculateMaf)

    FindHWEofSNPs = findHWEofSNPs(CalculateSnpSkewStatus)

    GenerateHwePlot  = generateHwePlot(FindHWEofSNPs)

    plots = channel.empty().mix(
        DifferentialMissingnessPlot,
        GenerateHwePlot,
        ShowHWEStats,
        ShowInitMAF,
        CalculateMaf,
        GenerateSnpMissingnessPlot,
        GenerateMissHetPlot).collect()

    collectPlotsTogetherAndZip(plots)
}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}