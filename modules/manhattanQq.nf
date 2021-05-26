def getInputChannels() {

	return channel.fromPath([
		"input/${params.cohortName}.bed",
		"input/${params.cohortName}.bim",
		"input/${params.cohortName}.fam"])
}

process getAssociationReport {
	container "quay.io/h3abionet_org/py3plink"
	input:
		path inputFiles
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

    subject = "[nextflow|h3agwaws] run ${workflow.runName} has finished"
    attachment1 = "${params.outputDir}/manhattanPlot.pdf"
    attachment2 = "${params.outputDir}/qqplot.pdf"
    message = \
        """\
        Hi there, 

        Your nextflow job ${workflow.scriptName}: ${workflow.runName} has finished.
        Please check the attachments to this email,
        and the execution summary below. 

        All the best,
        H 3 A G W A S



        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
    .stripIndent()

    sendMail(
        to: "${params.email}",
        subject: "${subject}",
        body: "${message}",
        attach: "${attachment1}", "${attachment2}")
}