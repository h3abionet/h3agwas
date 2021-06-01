include {
    userEmailAddressIsProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

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
            --bfile ${params.cohortName} \
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