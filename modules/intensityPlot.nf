def getInputChannels() {
    return channel.fromPath(
            params.inputFileDir +
            params.inputFilePrefix +
            '_gtReport_File-*')
}


process filterRecordsForChosenSnv() {
    input:
        path genotypeDataFile
    output:
        path "${genotypeDataFile}.xy.csv"
    script:
        """
        zcat ${genotypeDataFile} | grep "${params.snvName}," > ${genotypeDataFile}.xy.csv
        """
}

process concatenateDataFiles {
    input:
        path dataFiles
    output:
        path "${params.snvName}.csv"
    script:
        """
        echo 'SNP Name,Sample ID,Allele1 - Top,Allele2 - Top,GC Score,X,Y,B Allele Freq,Log R Ratio' > ${params.snvName}.csv
        cat ${params.inputFilePrefix}* >> ${params.snvName}.csv
        """
}

process plotXYintensityData {
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path genotypeDataForChosenSnv
    output:
        path "${params.snvName}_XYplot.pdf"
    script:
        """
        #!/usr/bin/env Rscript --vanilla
        library(tidyverse)
        snvdata <- read.csv(file="${genotypeDataForChosenSnv}")
        ggplot(snvdata, aes(x=X,y=Y)) +
            geom_point() +
            ggtitle("snv: ${params.snvName}")
        ggsave("${params.snvName}_XYplot.pdf")
        """
}

