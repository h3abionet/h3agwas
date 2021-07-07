include {
    checkOutputDir;
    checkCohortName;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    checkInputCohortData;
    getCohortData;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkInputCohortData('basicFiltered')
}

def getInputChannels() {
    return getCohortData('basicFiltered')
}

process getPopulationStratificationReports {
    label 'plink'
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        tuple path("${params.cohortName}.eigenval"), path("${params.cohortName}.eigenvec")
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --pca \
            --out ${params.cohortName}
        """
}

process drawPopulationStratificationPlot {
    label 'matplotlib'
    label 'mediumMemory'

    tag "eignvals, eigenvecs"

    input:
        tuple path(cohortEigenval), path(cohortEigenvec)
    output:
        publishDir "${params.outputDir}/sampleFiltered/plots", mode: 'copy'
        path "${params.cohortName}.pca.pdf"
        path "${params.cohortName}.eigenvalue.pdf"
    script:
        cc_fname = 0
        cc       = 0
        col      = 0
        output="${params.cohortName}.pca.pdf"
        template "drawPopulationStratificationPlot.py"
}

process getDiscordantSampleSexInfoReport {
    label 'plink'
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        publishDir "${params.outputDir}/sampleFiltered/reports", mode: 'copy'
        path "${params.cohortName}.sexcheck"

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --check-sex \
                ${params.sampleQC.maxInbreedingCoefficientForFemaleCalls} \
                ${params.sampleQC.minInbreedingCoefficientForMaleCalls} \
            --out ${params.cohortName}
        """
}

process selectSamplesWithDiscordantSexInfo {
    label 'smallMemory'

    tag "cohortData"

    input:
        path cohortSexcheck
    output:
        path "samplesWithDiscordantSexInfo.txt"
    script:
        """
        awk \
            '/PROBLEM/ {print \$1 "\t" \$2}' \
            ${cohortSexcheck} \
            > samplesWithDiscordantSexInfo.txt
        """
}

process getSampleBasedMissingnessReport {
    label 'plink'
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        publishDir "${params.outputDir}/sampleFiltered/reports", mode: 'copy'
        path "${params.cohortName}.imiss"

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --missing \
            --out ${params.cohortName}
        """
}

process drawSampleMissingnessHistogram {
    label 'matplotlib'
    label 'mediumMemory'

    tag "cohortImiss"

    input:
        path cohortImiss

    output:
        publishDir "${params.outputDir}/sampleFiltered/plots", mode: 'copy'
        path output

    script:
        input  = cohortImiss
        base   = cohortImiss.getBaseName()
        label  = "samples"
        output = "${params.cohortName}.sampleMissingness.pdf"
        template "drawSampleMissingnessHistogram.py"
}

process getIdentityByDescentReport {
    label 'plink'
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        publishDir "${params.outputDir}/sampleFiltered/reports", mode: 'copy'
        path "${params.cohortName}.genome"

    script:
        """
        plink \
            --keep-allele-order \
            -bfile ${cohortBed.getBaseName()} \
            --threads ${task.cpus} \
            --min ${params.sampleQC.minRelatednessPiHat} \
            --genome \
            --out ${params.cohortName}
        """
}

process selectSamplesWithHighRelatedness {
    label 'python3'
    label 'mediumMemory'

    tag "genome, imiss"

    input:
        path cohortGenome
        path cohortImiss

    output:
        path outfname

    script:
        base = cohortImiss.baseName
        outfname = "samplesWithHighRelatedness.txt"
        template "selectSamplesWithHighRelatedness.py"
}

process getSampleHeterozygosityReport {
    label 'plink'
    label 'mediumMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        publishDir "${params.outputDir}/sampleFiltered/reports", mode: 'copy'
        path "${params.cohortName}.het"
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --het \
            --out ${params.cohortName}
        """
}

process drawMissingnessHeterozygosityPlot {
    label 'matplotlib'
    label 'smallMemory'

    tag "imiss, het"

    input:
        path cohortImiss
        path cohortHet
    output:
        publishDir "${params.outputDir}/sampleFiltered/plots", mode: 'copy'
        path output
    script:
        base = cohortImiss.baseName
        output  = "${params.cohortName}.missingnessHeterozygosity.pdf"
        template "drawMissingnessHeterozygosityPlot.py"
}

process selectSamplesFailingMHTest {
    label 'python3'
    label 'mediumMemory'

    tag "imiss, het"

    input:
        path cohortImiss
        path cohortHet
    output:
        path outfname
    script:
        base = cohortImiss.baseName
        outfname = "samplesFailingMHTest.txt"
        template "selectSamplesFailingMHTest.py"
}

process concatenateSampleLists {

    tag "all sample lists"

    input:
        path sampleLists
    output:
        path "concatenatedSampleList.txt"
    script:
        """
        cat *.txt | sort -k2 | uniq > concatenatedSampleList.txt
        """
}

process removeLowQualitySamples {
    label 'plink'
    label 'smallMemory'

    tag "cohortData, lowQualitySamples"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
        path lowQualitySamples

    output:
        publishDir path: "${params.outputDir}/sampleFiltered/cohortData", mode: 'copy'
        path "${params.cohortName}.sampleFiltered.{bed,bim,fam}"

    script:
        """
         plink \
            --keep-allele-order \
            -bfile ${cohortBed.getBaseName()} \
            --remove ${lowQualitySamples} \
            --make-bed \
            --out ${params.cohortName}.sampleFiltered
        """
}

def sendWorkflowExitEmail() {
    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: getBasicEmailSubject(),
            body: getBasicEmailMessage(),
            attach: "${params.outputDir}plotArchives/sampleFiltering.tar.gz")
    }
}
