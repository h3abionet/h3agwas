include {
    userEmailAddressIsProvided;
    checkInputDir;
    checkCohortName;
    checkGenotypeReportPrefix;
    checkEmailAdressProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
    getCohortData;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkInputDir()
    checkCohortName()
    checkEmailAdressProvided()
}

def getInputChannels() {
	return getCohortData(params.inputStep)
}

process getAssociationReport {
	label 'plink'
	input:
		tuple path(cohortBed), path(cohortBim), path(cohortFam)
	output:
		path "${params.cohortName}.report.assoc"
	script:
		"""
		 plink \
		 	--bfile ${params.cohortName} \
		 	--assoc \
		 	--maf 0.05 \
		 	--out ${params.cohortName}.report 
		"""
}

process drawManhattanPlot {
	label 'qqman'
	publishDir "${params.outputDir}", mode: 'copy'
	input:
		path associationReport
	output:
		path manhattanPlot
	script:
		manhattanPlot = "manhattan.pdf"
		"""
		#!/usr/bin/env Rscript --vanilla
		library(qqman)
		assoc <- read.table("${associationReport}", header=TRUE)
		pdf("${manhattanPlot}")
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
	publishDir "${params.outputDir}", mode: 'copy'
	input:
		path associationReport
	output:
		path qqplot
	script:
		qqplot = "qqplot.pdf"
		"""
		#!/usr/bin/env Rscript --vanilla
		library(qqman)
		assoc <- read.table("${associationReport}", header=TRUE)
		pdf("${qqplot}")
		qq(assoc\$P)
		dev.off()
		"""
}

def sendWorkflowExitEmail() {

    subject = getBasicEmailSubject()
    attachment = [
        "${params.outputDir}manhattan.pdf",
        "${params.outputDir}qqplot.pdf"]
    message = getBasicEmailMessage()

    if (userEmailAddressIsProvided()) {
	    sendMail(
	        to: "${params.email}",
	        subject: "${subject}",
	        body: "${message}",
	        attach: ["${attachment[0]}", "${attachment[1]}"])
	}
}
