include {
    getCohortData;
    checkCohortName;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkEmailAdressProvided()
}

def getInputChannels() {
    return getCohortData('sampleFiltered')
}

process getMissingnessReport {
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "plink.missing"

    script:
        """
        plink \
            --bfile ${cohortBed.getBaseName()} \
            --threads ${task.cpus} \
            --autosome \
            --test-missing \
                mperm=${params.differentialMissingness.mperm}
        """
}

process drawDifferentialMissingnessPlot {
    label 'smallMemory'

    input:
        path cohortMissing

    output:
        path output

    script:
        input = cohortMissing
        output = "differentialMissingness.pdf"
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
            '{ if (\$5 <= ${params.differentialMissingness.cut}) { print \$2 } }' \
            ${cohortMissing} \
            > highDifferentialMissingnessSnvs.txt
        """
}

process getAlleleFrequencyReport {
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "plink.frq"

    script:
        """
        plink \
            --bfile ${cohortBed.getBaseName()} \
            --freq \
        """
}

process drawAlleleFrequencyPlot {
    label 'smallMemory'

    input:
        path cohortFrq

    output:
        path output

    script:
        input = cohortFrq
        output  = "alleleFrequencies.pdf"
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
            '\$5<${params.minMAF} { print \$2 }' \
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
            'NR!=1 && (\$3>=${params.missingness.maxRate.cases} || \$4>=${params.missingness.maxRate.controls}) { print \$2 }' \
            ${cohortMissing} \
            > highMissingnessSnvs.txt
        """
}

process getHardyWeinbergEquilibriumReport {
    label 'smallMemory'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "plink.hwe"

    script:
        """
        plink \
            --bfile ${cohortBed.getBaseName()} \
            --hardy
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
    label 'smallMemory'

    input:
        path controlOnlyCohortHwe

    output:
        path output

    script:
        input  = controlOnlyCohortHwe
        output = "controlOnlyCohortHwe.pdf"
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
            '{ if (\$9 <= ${params.hwe.cut}) { print \$2 } }' \
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

    publishDir path: "${params.outputDir}", mode: 'copy'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
        path lowQualitySnvs

    output:
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
          body: getBasicEmailMessage())
  }
}
