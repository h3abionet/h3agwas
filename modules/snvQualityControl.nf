include {
    getCohortData;
    checkCohortName;
    checkOutputDir;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
    checkInputCohortData;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkInputCohortData('sampleFiltered')
}

def getInputChannels() {
    return getCohortData('sampleFiltered')
}

process getMissingnessReport {
    label 'plink'
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        publishDir "${params.outputDir}/snvFiltered/reports", mode: 'copy'
        path "${params.cohortName}.missing"

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --threads ${task.cpus} \
            --autosome \
            --test-missing \
                mperm=${params.snvQC.numberOfMaxTPermutations} \
            --out ${params.cohortName}
        """
}

process drawDifferentialMissingnessPlot {
    label 'matplotlib'
    label 'mediumMemory'

    input:
        path cohortMissing

    output:
        publishDir "${params.outputDir}/snvFiltered/plots", mode: 'copy'
        path output

    script:
        input = cohortMissing
        output = "${params.cohortName}.differentialMissingness.pdf"
        template "drawDifferentialMissingnessPlot.py"
}

process selectSnvsWithHighDifferentialMissingness {
    label 'smallMemory'

    input:
        path cohortMissing

    output:
        path "highDifferentialMissingnessSnvs.txt"

    script:
        """
        awk \
            '{ if (\$5 <= ${params.snvQC.minDifferentialMissingnessPvalue}) { print \$2 } }' \
            ${cohortMissing} \
            > highDifferentialMissingnessSnvs.txt
        """
}

process getAlleleFrequencyReport {
    label 'plink'
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        publishDir "${params.outputDir}/snvFiltered/reports", mode: 'copy'
        path "${params.cohortName}.frq"

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --freq \
            --out ${params.cohortName}
        """
}

process drawAlleleFrequencyPlot {
    label 'matplotlib'
    label 'mediumMemory'

    input:
        path cohortFrq

    output:
        publishDir "${params.outputDir}/snvFiltered/plots", mode: 'copy'
        path output

    script:
        input = cohortFrq
        output  = "${params.cohortName}.alleleFrequencies.pdf"
        template "drawAlleleFrequencyPlot.py"
}

process selectlowMinorAlleleFrequencySnvs {
    label 'smallMemory'

    input:
        path cohortFrq
    output:
        path "lowMinorAlleleFrequencySnvs.txt"
    script:
        """
        awk \
            '\$5<${params.snvQC.minAlleleFrequency} { print \$2 }' \
            ${cohortFrq} \
            > lowMinorAlleleFrequencySnvs.txt
        """
}

process selectSnvsWithHighMissingness {
    label 'smallMemory'

    input:
        path cohortMissing

    output:
        path "highMissingnessSnvs.txt"

    script:
        """
        awk \
            'NR!=1 && (\$3>=${params.snvQC.maxMissingnessPerSnv.cases} || \$4>=${params.snvQC.maxMissingnessPerSnv.controls}) { print \$2 }' \
            ${cohortMissing} \
            > highMissingnessSnvs.txt
        """
}

process getHardyWeinbergEquilibriumReport {
    label 'plink'
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        publishDir "${params.outputDir}/snvFiltered/reports", mode: 'copy'
        path "${params.cohortName}.hwe"

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --hardy \
            --out ${params.cohortName}
        """
}

process selectControlOnlyHWETests {
    label 'smallMemory'

    tag 'cohortHwe'

    input:
        path cohortHwe

    output:
        path output

    script:
        output = "controlOnly.hwe"
        """
        head -1 ${cohortHwe} > $output
        grep 'UNAFF' ${cohortHwe} >> $output
        """
}

process drawHardyWeinbergEquilibriumPlot {
    label 'matplotlib'
    label 'mediumMemory'

    input:
        path controlOnlyCohortHwe

    output:
        publishDir "${params.outputDir}/snvFiltered/plots", mode: 'copy'
        path output

    script:
        input  = controlOnlyCohortHwe
        output = "${params.cohortName}.controlOnlyCohortHwe.pdf"
        template "drawHardyWeinbergEquilibriumPlot.py"
}

process selectSnvsOutOfHardyWeinbergEquilibrium {
    label 'smallMemory'

    input:
        path controlOnlyCohortHwe

    output:
        path "outOfHardyWeinbergEquilibriumSnvs.txt"

    script:
        """
        awk \
            '{ if (\$9 <= ${params.snvQC.minHardyWeinbergEquilibriumPvalue}) { print \$2 } }' \
            ${controlOnlyCohortHwe} \
            > outOfHardyWeinbergEquilibriumSnvs.txt        
        """
}

process concatenateSnvLists {

    tag "all snv lists"

    input:
        path sampleLists
    output:
        path "concatenatedSnvList.txt"
    script:
        """
        cat *.txt | sort -k1 | uniq > concatenatedSnvList.txt
        """
}

process removeLowQualitySnvs {
    label 'plink'
    label 'smallMemory'

    tag "cohortData, lowQualitySnvs"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
        path lowQualitySnvs

    output:
        publishDir path: "${params.outputDir}/snvFiltered/cohortData", mode: 'copy'
        path "${params.cohortName}.snvFiltered.{bed,bim,fam}"

    script:
        """
         plink \
            --keep-allele-order \
            -bfile ${cohortBed.getBaseName()} \
            --exclude ${lowQualitySnvs} \
            --make-bed \
            --out ${params.cohortName}.snvFiltered
        """
}

def sendWorkflowExitEmail() {
    if (userEmailAddressIsProvided()) {
      sendMail(
          to: "${params.email}",
          subject: getBasicEmailSubject(),
          body: getBasicEmailMessage(),
	  attach: "${params.outputDir}plotArchives/snvFiltering.tar.gz")
  }
}
