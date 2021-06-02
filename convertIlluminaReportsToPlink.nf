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
} from "${projectDir}/modules/intensityPlot.nf"

include {
    checkInputParams;
    getInputChannels;
    convertGenotypeReportsToLgen;
    concatenateLgenFiles;
    getFamFileFromSampleReport;
    getMapFileFromSnpReport;
    convertPlinkLongFormatToPlinkBinary
    sendWorkflowExitEmail;
} from "${projectDir}/modules/illuminaReportsToPlink.nf"


workflow {

    checkInputParams()

    (genotypeReports,
    sampleReport,
    snpReport) = getInputChannels()

    lgenFiles = convertGenotypeReportsToLgen( genotypeReports )
    lgenFile = concatenateLgenFiles( lgenFiles )
    famFile = getFamFileFromSampleReport( sampleReport )
    mapFile = getMapFileFromSnpReport( snpReport )

    convertPlinkLongFormatToPlinkBinary(
        lgenFile, 
        famFile, 
        mapFile)

}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}