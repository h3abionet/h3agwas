def checkOutputDir() {
    if (stringIsNull(params.outputDir)) {
        exit 1, 'params.outputDir not set -> please provide an output directory with enough space to save the results'
    }   
}
def checkSnvName() {
    if (stringIsNull(params.snvName)) {
        exit 1, 'please provide a snv id!'
    }
}
def checkIlluminaGenotypeReports() {
    if (stringIsNull(params.illumina.genotypeReports)) {
        exit 1, 'params.illumina.genotypeReports not set -> please provide a genotype report path glob pattern'
    }
    // also need to check if any file exists... //
}
def checkIlluminaSampleReport() {
    if (stringIsNull(params.illumina.sampleReport)) {
        exit 1, 'params.illumina.sampleReport not set -> please provide a sample report file path'
    }
    // also need to check if any file exists... //
}
def checkIlluminaLocusReport() {
    if (stringIsNull(params.illumina.locusReport)) {
        exit 1, 'params.illumina.locusReport not set -> please provide a locus report file path'
    }
    // also need to check if any file exists... //
}
def checkClinicalPhenotypeFam() {
    if (stringIsNull(params.clinicalPhenotypeFam)) {
        exit 1, 'params.clinicalPhenotypeFam not set -> please provide a clinical phenotype fam file path'
    }
    // also need to check if any file exists... //
}
def checkEmailAdressProvided() {
    if (!userEmailAddressIsProvided()) {
        println 'You have not specified an email address; ' \
        	+ 'we will not email you the results of this workflow.'
    }
}
def checkCohortName () {
	if (stringIsNull(params.cohortName)) {
		exit 1, 'params.cohortName not set -> please provide a short cohort name to label your output files'
	}
}
def checkReferencePanelsDir() {
    if (stringIsNull(params.referencePanelsDir)) {
        exit 1, 'please provide a directory of reference panels e.g. from 1000 Genomes'
    }
}
def checkReferenceSequence() {
    if (stringIsNull(params.referenceSequence)) {
        exit 1, 'params.referenceSequence not set -> please provide a reference sequence fasta file path'
    }
    // also need to check if any file exists... //
}

def userEmailAddressIsProvided() {
	return !(stringIsNull(params.email))
}
def stringIsNull(string) {
	return ( string =~ /NULL/ )
}

def getBasicEmailSubject() {
    return "[nextflow|h3agwaws] run ${workflow.runName} has finished"
}
def getBasicEmailMessage() {
    return """\
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
}

def getCohortData(inputDataTag) {

    dataTag = (inputDataTag == '') ? '' : ".${inputDataTag}"

    bed = channel.fromPath(
        params.outputDir + params.cohortName + "${dataTag}.bed")
    bim = channel.fromPath(
        params.outputDir + params.cohortName + "${dataTag}.bim")
    fam = channel.fromPath(
        params.outputDir + params.cohortName + "${dataTag}.fam")

    return bed.combine(bim).combine(fam)
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

process collectPlotsTogetherAndZip {

    input:
        val label
        path plots

    output:
        publishDir "${params.outputDir}", mode: 'copy'
        path "plots-${label}.tar.gz"

    script:
        """
        mkdir plots-${label}
        cp -L *.pdf plots-${label}
        tar -zcvf plots-${label}.tar.gz plots-${label}
        """
}