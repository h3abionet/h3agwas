#!/usr/bin/env nextflow

/*
 *  CONVERT GENOTYPE REPORTS TO PLINK BINARY {bed+bim+fam} FILESET
 *  ==============================================================
 *
 *  This script takes as input an illumina gsgt format dataset and  
 *  converts it to an lgen format of the PLINK long-format fileset
 *  and then to PLINK binary fileset: plink.bed, plink.bim, plink.fam
 *
 *  The gsgt illumina data format is a set of
 *  genotype reports. A genotpye report is an csv file that contains 
 *  a portion of the genotype results for the cohort. 
 *  There are several genotpye reports for a illumina sequencing 
 *  project, and one file does not correspond to any one snv or any 
 *  one sample.
 *
 *  The path to all the input genotype reports are specified in an 
 *  input tsv file, and this file is passed to the workflow by adding:
 *  --input <relative path to tsv file>
 *  to the command line.
 *
 ********************************************************************/

nextflow.enable.dsl = 2

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/base.nf"

include {
    checkInputParams;
    getInputChannels;
    getChunksFromGenotypeReports;
    concatenateLgenChunks;
    convertGenotypeReportsToLgen;
    getSnpReport;
    getSampleReport;
    getFamFileFromSampleReport;
    getMapFileFromSnpReport;
    convertPlinkLongFormatToPlinkBinary;
    alignGenotypeDataToReference;
    filterSitesWithoutRefOrAltAlleles;
    getFinalPlinkBinaryFileset;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/illuminaReportsToPlink.nf"

workflow {

    checkInputParams()

    (genotypeReports, 
     sampleReport, 
     locusReport, 
     phenotypeFam,
     referenceSequence) \
        = getInputChannels()

    genotypeReportChunks = splitTextFiles(genotypeReports)

    lgenChunks = convertGenotypeReportToLongFormat(genotypeReportChunks)

    cohortLgen = concatenateLgenFiles(lgenChunks)

    cohortMap = convertLocusReportToMap(locusReport)

    filteredSampleReport = removeSamplesWithFailedGenotypes(sampleReport)

    illuminaFam = convertSampleReportToFam(filteredSampleReport)

    cohortFam = intersectFamFiles(phenotypeFam, illuminaFam)

    cohortData = buildCohortData(cohortLgen, cohortMap, cohortFam)

    alignedGenotypeSet = alignGenotypesToReference(cohortData, referenceSequence)

    filteredGenotypeSet = selectBiallelicSnvs(alignedGenotypeSet)

    alignedCohortData = rebuildCohortData(filteredGenotypeSet, cohortFam)

}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}

