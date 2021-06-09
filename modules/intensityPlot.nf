include {
    userEmailAddressIsProvided;
    checkSnvName;
    checkInputDir;
    checkGenotypeReportPrefix;
    checkEmailAdressProvided;
    getBasicEmailMessage;
    getBasicEmailSubject;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkSnvName()
    checkInputDir()
    checkGenotypeReportPrefix()
    checkEmailAdressProvided()
}

def getInputChannels() {
    return channel.fromPath(
            params.inputDir +
            params.genotypeReportPrefix +
            '*')
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
            = "SNP Name,Sample ID,Allele1 - Top," \
            + "Allele2 - Top,GC Score,X,Y,B Allele Freq,Log R Ratio"
        """
        echo ${mergedReportHeader} > ${params.snvName}.csv
        cat ${params.genotypeReportPrefix}* >> ${params.snvName}.csv
        """
}

process drawXYintensityPlot {
    label "tidyverse"
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

def printWorkflowExitMessage() {
    if (workflow.success) {
        log.info "Workflow completed without errors".center(60)
    } else {
        log.error "Oops .. something went wrong!".center(60)
    }
    log.info "Check output files in folder:".center(60)
    log.info "${params.outputDir}".center(60)
}

def sendWorkflowExitEmail() {

    subject = getBasicEmailSubject()
    attachment = "${params.outputDir}/${params.snvName}_XYintensities.pdf"
    message = getBasicEmailMessage()

    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: "${subject}",
            body: "${message}",
            attach: "${attachment}")
    }
}
