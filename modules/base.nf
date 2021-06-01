def checkSnvName() {
    if (stringIsNull(params.snvName)) {
        exit 1, 'please provide a snv id!'
    }
}
def checkInputDir() {
    if (stringIsNull(params.inputDir)) {
        exit 1, 'please provide an input directory!'
    }
    // also need to check is input directory exists... //
}
def checkGenotypeReportPrefix() {
    if (stringIsNull(params.genotypeReportPrefix)) {
        exit 1, 'please provide a genotype report prefix!'
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
		exit 1, 'please provide a cohort name!'
	}
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