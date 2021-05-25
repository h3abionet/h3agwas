/*
 *	DRAW INTENSITY PLOT FOR A SNV
 *	=============================
 *
 *	This script takes as input an illumina gsgt format dataset, and a 
 *	SNV id chosen by the user, and then draws a plot of the intensity 
 *  (X,Y) values across all samples in the dataset for this chosen
 *  SNV.
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
nextflow.enable.dsl=2

include {
    getInputChannels;
    filterRecordsForChosenSnv;
    mergeGenotypeReports;
    plotXYintensityFields;
} from "${projectDir}/modules/intensityPlot.nf"

workflow {

    genotypeReports \
        = getInputChannels()

    genotypeReportsOfChosenSnv \
        = filterRecordsForChosenSnv(
            genotypeReports) 

    snvGenotypeReport \
        = mergeGenotypeReports(
            genotypeReportsOfChosenSnv.collect())

    plotXYintensityFields(
        snvGenotypeReport)

}
