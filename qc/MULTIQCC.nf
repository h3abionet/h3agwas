#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {

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

} from './modules/MULTIQCC.nf';

bed = channel.fromPath(
    "${params.input_dir}.bed")
bim = channel.fromPath(
    "${params.input_dir}.bim")
fam = channel.fromPath(
    "${params.input_dir}.fam")
imiss = channel.fromPath(
    "${params.input_dir}.imiss")
poor = channel.fromPath(
    "${params.input_dir}.txt")
plates = channel.fromPath(
    "${params.inputDir}plates")
poorgc10 = channel.fromPath(
    "${params.inputDir}poorgc10.lst")

Channel.fromPath(params.case_control).set { cc_ch }

// bed.combine(bim).combine(fam)                     

workflow {

    DuplicatePositions = getListOfDuplicatePositions(bed.combine(bim).combine(fam))

    RemovePoorSamples = removeSamplesWithPoorClinicalData(bed.combine(bim).combine(fam).combine(poor))

    RemoveDuplicateSnps = removeDuplicatedSnvPositions(RemovePoorSamples,DuplicatePositions)

    IdentifyIndivDiscSexinfo = identifyIndivDiscSexinfo(RemoveDuplicateSnps)

    GenerateSnpMissingnessPlot = generateSnpMissingnessPlot(IdentifyIndivDiscSexinfo)

    GenerateIndivMissingnessPlot = generateIndivMissingnessPlot(IdentifyIndivDiscSexinfo)

    GetInitMAF = getInitMAF(bed.combine(bim).combine(fam))

    ShowInitMAF = showInitMAF(GetInitMAF)

    ShowHWEStats = showHWEStats(IdentifyIndivDiscSexinfo)

    RemoveQCPhase1 = removeQCPhase1(bed.combine(bim).combine(fam))

    PruneForIBD  = pruneForIBD(RemoveQCPhase1)

    FindRelatedIndiv = findRelatedIndiv(IdentifyIndivDiscSexinfo,PruneForIBD)

    SampleHeterozygosity = calculateSampleHeterozygosity(RemoveQCPhase1)

    GenerateMissHetPlot = generateMissHetPlot(SampleHeterozygosity)

    BadIndivsMissingHet = getBadIndivsMissingHet(SampleHeterozygosity)

    NoSampleSheet = noSampleSheet()

    RemoveQCIndivs = removeQCIndivs(BadIndivsMissingHet,FindRelatedIndiv,IdentifyIndivDiscSexinfo,
                                    NoSampleSheet, RemoveQCPhase1)

    CalculateSnpSkewStatus = calculateSnpSkewStatus(RemoveQCIndivs.combine(cc_ch))

    DifferentialMissingnessPlot = generateDifferentialMissingnessPlot(CalculateSnpSkewStatus)

    FindSnpExtremeDifferentialMissingness = findSnpExtremeDifferentialMissingness(CalculateSnpSkewStatus)

    RemoveSkewSnps = removeSkewSnps(RemoveQCIndivs,FindSnpExtremeDifferentialMissingness)

    CalculateMaf = calculateMaf(RemoveSkewSnps)

    GenerateMafPlot = generateMafPlot(CalculateMaf)

    FindHWEofSNPs = findHWEofSNPs(CalculateSnpSkewStatus)

    GenerateHwePlot  = generateHwePlot(FindHWEofSNPs)

}