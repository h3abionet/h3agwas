include {
    checkOutputDir;
    checkCohortName;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getCohortData;
    getBasicEmailSubject;
    getBasicEmailMessage;
    checkReferenceSequence;
    checkInputCohortData;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkOutputDir()
    checkEmailAdressProvided()
    checkInputCohortData('input')
    checkReferenceSequence()
}

def getInputChannels() {
    return [
        getCohortData('input'),
        getReferenceSequence()]
}

def getReferenceSequence() {
    return channel
        .fromPath(params.baseQC.referenceSequence)
}

process selectDuplicatedVariants {
    label 'plink'

    tag "cohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path 'plink.dupvar'

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            --list-duplicate-vars \
            ids-only \
            suppress-first
        """
}

process removeDuplicatedVariants {
    label 'plink'

    tag "cohortData, duplicatedVariants"

    input:
    tuple path(cohortBed), path(cohortBim), path(cohortFam)
    path duplicatedVariants

    output:
        path("duplicatedVariantsRemoved.{bed,bim,fam}")

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()} \
            -exclude ${duplicatedVariants} \
            --make-bed \
            --out duplicatedVariantsRemoved
        """
}


process alignGenotypesToReference() {
    label 'mediumMemory'
    label 'plink2'

    tag "filteredCohortData, reference"

    cache 'lenient'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
        path referenceSequence
    output:
        path "temporary.vcf.gz"
    script:
        plinkBase = cohortBed.getBaseName()
        """
        plink2 \
            --bfile ${plinkBase} \
            --fa ${referenceSequence} \
            --ref-from-fa force \
            --normalize \
            --threads $task.cpus \
            --export vcf-4.2 id-paste=iid bgz \
            --real-ref-alleles \
            --out temporary
        """
}

process selectBiallelicSnvs() {
    label 'mediumMemory'
    label 'bcftools'

    tag "alignedGenotypes"

    input:
        path genotypeSet
    output:
        path "filtered.vcf.gz"
    script:
        """
        bcftools \
            view \
            -m2 \
            -M2 \
            -v snps \
            --threads $task.cpus \
            -Oz \
            -o filtered.vcf.gz \
        ${genotypeSet}
        """
}

process rebuildCohortData() {
    label 'mediumMemory'
    label 'plink2'

    tag "alignedGenotypes, filteredCohortFam"

    input:
        path alignedGenotypes
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
    output:
        path "aligned.{bed,bim,fam}"
    script:
        """
        plink2 \
            --vcf ${alignedGenotypes} \
            --fam ${cohortFam} \
            --threads $task.cpus \
            --make-bed \
            --double-id \
            --out aligned
        """
}

process removeReallyLowQualitySamplesAndSnvs {
    label 'plink'
    label 'smallMemory'

    tag "alignedCohortData"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        publishDir "${params.outputDir}/basicFiltered/cohortData", mode: 'copy'
        path "${params.cohortName}.basicFiltered.{bed,bim,fam}"

    script:
        """
        plink \
            --keep-allele-order \
            --bfile ${cohortBed.getBaseName()}  \
            --mind ${params.baseQC.maxMissingnessPerSample} \
            --make-bed \
            --out temp1
        plink \
            --keep-allele-order \
            --bfile temp1 \
            --geno ${params.baseQC.maxMissingnessPerSnv} \
            --make-bed \
            --out temp2
        plink \
            --keep-allele-order \
            --bfile temp2 \
            --maf ${params.baseQC.minAlleleFrequency} \
            --make-bed \
            --out temp3
        plink \
            --keep-allele-order \
            --bfile temp3  \
            --hwe ${params.baseQC.minHardyWeinbergEquilibriumPvalue} \
            --make-bed \
            --out ${params.cohortName}.basicFiltered     
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
