include {
    checkCohortName;
    userEmailAddressIsProvided;
    checkOutputDir;
    checkSelectedSnv;
    checkEmailAdressProvided;
    checkIlluminaGenotypeReports;
    getBasicEmailMessage;
    getBasicEmailSubject;
} from "${projectDir}/modules/base.nf"

include {
    getIlluminaGenotypeReports;
} from "${projectDir}/modules/buildInput.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkIlluminaGenotypeReports()
    checkSelectedSnv()
}

def getInputChannels() {
    return getIlluminaGenotypeReports()
}

process filterRecordsForChosenSnv() {
    input:
        path genotypeReport
    output:
        path "${genotypeReport}.xy.csv"
    script:
        """
        zcat ${genotypeReport} \
            | grep "${params.selectedSnv}," \
            > ${genotypeReport}.xy.csv
        """
}

process mergeGenotypeReports {
    input:
        path genotypeReportList
    output:
        path "${params.selectedSnv}.csv"
    script:
        mergedReportHeader \
            = "SNP Name,Sample ID,Allele1 - Top," \
            + "Allele2 - Top,GC Score,X,Y,B Allele Freq,Log R Ratio"
        """
        echo ${mergedReportHeader} > ${params.snvName}.csv
        cat ${params.genotypeReportPrefix}* >> ${params.selectedSnv}.csv
        """
}

process drawXYintensityPlot {
    label "tidyverse"

    input:
        path snvGenotypeReport
    output:
        publishDir "${params.outputDir}/intensityPlots", mode: 'copy'
        path "${params.cohortName}.${params.selectedSnv}.XY.pdf"
    script:
        """
        #!/usr/bin/env Rscript --vanilla
        library(tidyverse)
        snvGenotypeReport <- read.csv(file="${snvGenotypeReport}")
        ggplot(snvGenotypeReport, aes(x=X,y=Y)) +
            geom_point() +
            ggtitle("snv: ${params.selectedSnv}")
        ggsave("${params.cohortName}.${params.selectedSnv}.XY.pdf")
        """
}

def sendWorkflowExitEmail() {
    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: getBasicEmailSubject(),
            body: getBasicEmailMessage(),
            attach: "${params.outputDir}/intensityPlots/${params.cohortName}.${params.selectedSnv}.XY.pdf")
    }
}
