#!/usr/bin/env nextflow

/*
 *  CONVERT GENOTYPE REPORTS TO LGEN
 *  ================================
 *
 *  This script takes as input an illumina gsgt format dataset and  
 *  converts it to an lgen format of the PLINK long-format fileset
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

   getGenotypeReports;
   
   getSampleReport;
   
   getSnpReport;
   
   convertGenotypeReportsToLgen;
   
   concatenateLgenFiles;
   
   getFamFileFromSampleReport;
   
   getMapFileFromSnpReport;
   
   convertPlinkLongFormatToPlinkBinary

} from "${projectDir}/modules/illuminaReportsToPlink.nf"


workflow {

   genotypeReports = getGenotypeReports()

   sampleReport = getSampleReport()
   
   snpReport = getSnpReport()

   
   famFile = getFamFileFromSampleReport( sampleReport )
   
   mapFile = getMapFileFromSnpReport( snpReport )
   
   lgenFiles = convertGenotypeReportsToLgen( genotypeReports )
   
   lgenFile = concatenateLgenFiles( lgenFiles )
   
   convertPlinkLongFormatToPlinkBinary( lgenFile, 
                                        famFile, 
                                        mapFile )

}
