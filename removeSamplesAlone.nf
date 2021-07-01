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
    generateIndivMissingnessPlot;
    removeQCPhase1;
    pruneForIBD;
    findRelatedIndiv;
    calculateSampleHeterozygosity;
    generateMissHetPlot;
    getBadIndivsMissingHet;
    noSampleSheet;
    removeQCIndivs;
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

    GenerateIndivMissingnessPlot = generateIndivMissingnessPlot(IdentifyIndivDiscSexinfo)

    RemoveQCPhase1 = removeQCPhase1(cohortData)

    PruneForIBD  = pruneForIBD(RemoveQCPhase1)

    FindRelatedIndiv = findRelatedIndiv(IdentifyIndivDiscSexinfo,PruneForIBD)

    SampleHeterozygosity = calculateSampleHeterozygosity(RemoveQCPhase1)

    GenerateMissHetPlot = generateMissHetPlot(SampleHeterozygosity)

    BadIndivsMissingHet = getBadIndivsMissingHet(SampleHeterozygosity)

    NoSampleSheet = noSampleSheet()

    RemoveQCIndivs = removeQCIndivs(BadIndivsMissingHet,FindRelatedIndiv,IdentifyIndivDiscSexinfo,
                                    NoSampleSheet, RemoveQCPhase1)

    plots = channel.empty().mix(
        GenerateIndivMissingnessPlot,	
        GenerateMissHetPlot).collect()

    collectPlotsTogetherAndZip(plots)
}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}
