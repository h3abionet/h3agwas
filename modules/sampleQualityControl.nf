include {
    checkOutputDir;
    checkCohortName;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getCohortData;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkOutputDir()
    checkCohortName()
    checkEmailAdressProvided()
}

def getInputChannels() {
    return getCohortData('input')
}

process getDiscordantSampleSexInfoReport {
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "plink.sexcheck"

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --check-sex \
                ${params.checkSex.femaleMaxF} \
                ${params.checkSex.maleMinF}
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
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "plink.imiss"

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --missing
        """
}

process drawSampleMissingnessHistogram {
    label 'smallMemory'

    tag "cohortImiss"

    input:
        path cohortImiss
  
    output:
        path output
  
    script:
        input  = cohortImiss
        base   = cohortImiss.getBaseName()
        label  = "samples"
        output = "${base}-indmiss_plot".replace(".","_")+".pdf"
        template "drawSampleMissingnessHistogram.py"
}

process getIdentityByDescentReport {
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "plink.genome"

    script:
        """
        plink \
            --keep-allele-order \
            -bfile ${cohortBed.getBaseName()} \
            --threads ${task.cpus} \
            --min ${params.relatedness.piHat} \
            --genome
        """
}

process selectSamplesWithHighRelatedness {
    label 'smallMemory'

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
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        path "plink.het"
    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --het
        """
}

process drawMissingnessHeterozygosityPlot {
    label 'plink'
    label 'smallMemory'

    tag "imiss, het"

    input:
        path cohortImiss
        path cohortHet
    output:
        path output
    script:
        base = cohortImiss.baseName
        output  = "${base}-imiss-vs-het".replace(".","_")+".pdf"
        template "drawMissingnessHeterozygosityPlot.py"
}

process selectSamplesFailingMHTest {
    label 'smallMemory'

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

    publishDir path: "${params.outputDir}", mode: 'copy'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
        path lowQualitySamples

    output:
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

process collectPlotsTogetherAndZip {

    input:
        path plots

    output:
        publishDir "${params.outputDir}", mode: 'copy'
        path "plots-sampleFiltering.tar.gz"

    script:
        """
        mkdir plots-sampleFiltering
        cp -L *.pdf plots-sampleFiltering
        tar -zcvf plots-sampleFiltering.tar.gz plots-sampleFiltering
        """
}


def sendWorkflowExitEmail() {
    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: getBasicEmailSubject(),
            body: getBasicEmailMessage(),
            attach: "${params.outputDir}/quality-control/sampleA-Samples-plots.tar.gz")
    }
}
