include {
    checkInputDir;
    checkCohortName;
    checksamplesWithPoorClinicalData;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkInputDir()
    checkCohortName()
    checkEmailAdressProvided()
    checksamplesWithPoorClinicalData()
}

def getInputChannels() {

    bed = channel.fromPath(
        "${params.inputDir}${params.cohortName}.bed")
    bim = channel.fromPath(
        "${params.inputDir}${params.cohortName}.bim")
    fam = channel.fromPath(
        "${params.inputDir}${params.cohortName}.fam")

    return [
        bed.combine(bim).combine(fam),
        channel.fromPath("${params.samplesWithPoorClinicalData}")]

}

process getListOfDuplicatePositions {
    label 'plink'
    input:
        tuple path(bed), path(bim), path(fam)
    output:
        tuple path('plink.dupvar'), path('plink.hh')
    script:
        """
        plink \
            --bfile ${params.cohortName} \
            --list-duplicate-vars \
            ids-only \
            suppress-first
        """
}

process removeSamplesWithPoorClinicalData {
    label 'plink'
    input:
        tuple path(bed), path(bim), path(fam)
        path(poorSamples)
    output:
        path("${output}.{bed,bim,fam}")        
    script:
        output = "${params.cohortName}.sampleFiltered"
        """
        plink \
            --bfile ${params.cohortName} \
            --remove ${poorSamples} \
            --make-bed \
            --out ${output}
        """
}

process removeDuplicatedSnvPositions {
    label 'plink'
    publishDir "${params.outputDir}/quality-control", mode: 'copy'
    input:
        tuple path(bed), path(bim), path(fam)
        tuple path('plink.dupvar'), path('plink.hh')
    output:
        path("${output}.{bed,bim,fam,log}")
    script:
        output = "${params.cohortName}.DuplicatesRemoved"
        """
        plink \
            --bfile ${params.cohortName}.sampleFiltered \
            --make-bed \
            -exclude plink.dupvar \
            --out ${output}
        """
}

def sendWorkflowExitEmail() {

    subject = getBasicEmailSubject()
    message = getBasicEmailMessage()

    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: "${subject}",
            body: "${message}")
    }
}