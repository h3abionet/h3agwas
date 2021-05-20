/*
 *	DRAW INTENSITY PLOT FOR A SNV
 *	=============================
 *
 *	This script will take as input an illumina gsgt dataset, and a 
 *	SNV id, and then draw a plot of the intensity (X,Y) values across
 *	all samples in the dataset for this fixed SNV. 
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    getInputChannels;
    filterRecordsForChosenSnv;
    concatenateDataFiles;
    plotXYintensityData;
} from "${projectDir}/modules/intensityPlot.nf"

workflow {

    genotypeDataFiles\
        = getInputChannels()

    genotypeDataFilesForChosenSnv\
        = filterRecordsForChosenSnv(
            genotypeDataFiles) 

    genotypeDataForChosenSnv\
        = concatenateDataFiles(
            genotypeDataFilesForChosenSnv.collect())

    plotXYintensityData(
        genotypeDataForChosenSnv)

}
