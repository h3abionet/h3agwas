include {
    userEmailAddressIsProvided;
    checkOutputDir;
    checkCohortName;
    checkEmailAdressProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
    getCohortData;
    checkInputCohortData;
    checkAssociationInput;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkAssociationInput()
    checkInputCohortData(params.associationInput)
}

def getInputChannels() {
    return getCohortData(params.associationInput)
}

process getAssociationReport {
    label 'plink'

    tag "${params.associationInput}CohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        publishDir "${params.outputDir}/${params.associationInput}/reports", mode: 'copy'
        path "${params.cohortName}.assoc"
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
	    --assoc \
	    --maf 0.01 \
	    --out ${params.cohortName}
	"""
}

process drawManhattanPlot {
    label 'qqman'

    tag "associationReport"

    input:
        path associationReport
    output:
        publishDir "${params.outputDir}/${params.associationInput}/plots", mode: 'copy'
        path manhattanPlot
    script:
        manhattanPlot = "${params.cohortName}.manhattan.png"
        """
        #!/usr/bin/env Rscript --vanilla
        library(qqman)
        assoc <- read.table("${associationReport}", header=TRUE)
        png("${manhattanPlot}", width = 1000, height = 1000, units = "px")
        manhattan(
	    assoc,
	    chr="CHR",
	    bp="BP",
	    snp="SNP",
	    p="P",
	    logp=TRUE,
	    ylim = c(0, 8))
        dev.off()
        """
}

process drawQqPlot {
    label 'qqman'
	
    tag "associationReport"

    input:
        path associationReport
    output:
        publishDir "${params.outputDir}/${params.associationInput}/plots", mode: 'copy'
        path qqplot
    script:
        qqplot = "${params.cohortName}.qqplot.png"
        """
        #!/usr/bin/env Rscript --vanilla
        library(qqman)
        assoc <- read.table("${associationReport}", header=TRUE)
        png("${qqplot}", width = 1000, height = 1000, units = "px")
        qq(assoc\$P)
        dev.off()
        """
}

def sendWorkflowExitEmail() {
    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: getBasicEmailSubject(),
            body: getBasicEmailMessage(),
            attach: "${params.outputDir}/plotArchives/association-${params.associationInput}.tar.gz")
    }
}
