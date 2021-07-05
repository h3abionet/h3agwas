#!/usr/bin/env nextflow
/*
 *  BUILD INPUT COHORT DATA FOR ANALYSIS
 *  ====================================
 *
 *  This script builds the plink binary data files - bed, bim, fam - 
 *  needed in downstream analysis from illumina genotype reports and
 *  a clinical phenotype fam file. We refer to the plink binary files
 *  collectively as the *cohortData*. 
 *
 *  The input illumina genotype files required are:
 *      + genotype reports: n files containing all cohort genotypes
 *      + sample report: a file with details of the samples
 *      + locus report: a file with details of the loci
 * 
 *  The clinical phenotype file must be in *fam* format, and you
 *  must prepare it yourself before you start the analysis.
 *  It is neccesary that the sample IDs in the clinical fam file match
 *  those in the illumina sample report, *for the samples that are 
 *  present in both files*.
 *
 *  We start by pulling in the illumina reports, and the clinical
 *  phenotype fam file. The illumina locus report is converted to a 
 *  map file; the illumina sample report is converted to a fam file
 *  with blank phenotype information. We compute the intersection
 *  between the illumina fam and the clinical phenotype fam files; only
 *  sample records whose sample id is in both fam files are kept for
 *  downstream analysis, and the sample records *kept* are those from the 
 *  clinical phenotype fam, not the illumina fam. This means that both
 *  the self-reported sex and the phenotype values carried forward are
 *  those from the clinical fam; those values from the illumina fam are
 *  discarded.
 *
 *  We split the genotype reports into small chunks, and convert each
 *  chunk into a plink long-format genotype (lgen) file. These files
 *  are concatenated into a single lgen file for the cohort. We then
 *  use plink to build the cohortData from the map, fam, and lgen 
 *  files efficiently. 
 *
 ********************************************************************/
nextflow.enable.dsl = 2

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/base.nf"

include {
    checkInputParams;
    getInputChannels;
    splitTextFiles;
    convertGenotypeReportToLongFormat;
    concatenateLgenFiles;
    convertLocusReportToMap;
    removeSamplesWithFailedGenotypes;
    convertSampleReportToFam;
    intersectFamFilesBySampleId;
    buildCohortData;
    rebuildCohortDataWithPhenotypes;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/buildInput.nf"

workflow {

    checkInputParams()

    (genotypeReports,
     sampleReport,
     locusReport,
     phenotypeFam) \
        = getInputChannels()

    genotypeReportChunks \
        = splitTextFiles(
            genotypeReports)

    lgenChunks \
        = convertGenotypeReportToLongFormat(
            genotypeReportChunks)

    cohortLgen \
        = concatenateLgenFiles(
            lgenChunks)

    cohortMap \
        = convertLocusReportToMap(
            locusReport)

    filteredSampleReport \
        = removeSamplesWithFailedGenotypes(
            sampleReport)

    illuminaFam \
        = convertSampleReportToFam(
            filteredSampleReport)

    cohortData \
        = buildCohortData(
            cohortLgen,
            cohortMap,
            illuminaFam)

    cohortFam \
        = intersectFamFilesBySampleId(
            phenotypeFam,
            illuminaFam)

    phenotypedCohortData \
        = rebuildCohortDataWithPhenotypes(
            cohortData,
            cohortFam)

}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}

