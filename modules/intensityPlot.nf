def getInputChannels() {
    return channel.fromPath(
            params.inputFileDir +
            params.inputFilePrefix +
            '_gtReport_File-*')
}


process filterRecordsForChosenSnv() {
    input:
        path genotypeReport
    output:
        path "${genotypeReport}.xy.csv"
    script:
        """
        zcat ${genotypeReport} \
            | grep "${params.snvName}," \
            > ${genotypeReport}.xy.csv
        """
}

process mergeGenotypeReports {
    input:
        path genotypeReportList
    output:
        path "${params.snvName}.csv"
    script:
        mergedReportHeader \
            = "SNP Name,Sample ID,Allele1 - Top,"
            + "Allele2 - Top,GC Score,X,Y,B Allele Freq,Log R Ratio"
        """
        echo ${mergedReportHeader} > ${params.snvName}.csv
        cat ${params.inputFilePrefix}* >> ${params.snvName}.csv
        """
}

process plotXYintensityFields {
    publishDir "${params.outputDir}", mode: 'copy'

    input:
        path snvGenotypeReport
    output:
        path "${params.snvName}_XYintensities.pdf"
    script:
        """
        #!/usr/bin/env Rscript --vanilla
        library(tidyverse)
        snvGenotypeReport <- read.csv(file="${snvGenotypeReport}")
        ggplot(snvGenotypeReport, aes(x=X,y=Y)) +
            geom_point() +
            ggtitle("snv: ${params.snvName}")
        ggsave("${params.snvName}_XYintensities.pdf")
        """
}

