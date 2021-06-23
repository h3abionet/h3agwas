/*
 *  PHASE COHORT GENOTYPE DATA
 *  ==========================
 *
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    checkCohortName;
    checkReferencePanelsDir;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

workflow {

    checkInputParams()

    (cohortData, referencePanels) = getInputChannels()

    filteredReferencePanels = keepOnlyBiallelicSnvs(referencePanels)

    cohortGenotypes = convertPlinkBinaryToVcf(cohortData)

    alignedCohortGenotypesPerChromosome \
        = alignGenotypesToReference(
            cohortGenotypes.combine(filteredReferencePanels))

    indexedGenotypesPerChromosome = indexGenotypes(alignedCohortGenotypesPerChromosome)

    alignedCohortGenotypes \
        = concatenateAndSort(
            indexedGenotypesPerChromosome.collect())
}

def checkInputParams() {

    checkCohortName()
    checkReferencePanelsDir()
}

def getInputChannels() {

    inputDataTag = 'clean'

    bed = channel.fromPath(
        params.outputDir + params.cohortName + "-${inputDataTag}.bed")
    bim = channel.fromPath(
        params.outputDir + params.cohortName + "-${inputDataTag}.bim")
    fam = channel.fromPath(
        params.outputDir + params.cohortName + "-${inputDataTag}.fam")

    referencePanels = channel
        .of(1..22,'X')
        .map{ it -> [
            it,
            file(params.referencePanelsDir + '*chr' + it + '.*.vcf.gz')[0],
            file(params.referencePanelsDir + '*chr' + it + '.*.vcf.gz.tbi')[0]]}

    return [
        bed.combine(bim).combine(fam),
        referencePanels]
}

process keepOnlyBiallelicSnvs {
    label 'bcftools'
    label 'mediumMemory'

    tag "chromosome ${chromosome}"

    input:
        tuple val(chromosome), path(referencePanel), path(referencePanelIndex)
    output:
        tuple val("${chromosome}"), path("chr${chromosome}.refpanel.biallelic.vcf.gz")
    script:
        """
        bcftools \
            view \
            -M2 \
            -v snps \
            --threads $task.cpus \
            -Oz \
            -o chr${chromosome}.refpanel.biallelic.vcf.gz \
            ${referencePanel}
        """
}


process convertPlinkBinaryToVcf {
    label 'plink2'
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "${cohortBed.getBaseName()}.vcf.gz"

    script:
        """
        plink2 \
            --bfile ${cohortBed.getBaseName()} \
            --export vcf-4.2 bgz id-paste='fid' \
            --out ${cohortBed.getBaseName()}
        """
}

process alignGenotypesToReference {
    label 'beagle'
    label 'bigMemory'

    tag "chromosome ${chromosome}"

    input:
        tuple path(cohortGenotypes), val(chromosome), path(referencePanel)
    output:
        path "${cohortGenotypes.getSimpleName()}.${chromosome}.vcf.gz"
    script:
        """
        conform-gt \
            ref=${referencePanel} \
            gt=${cohortGenotypes} \
            chrom=${chromosome} \
            out="${cohortGenotypes.getSimpleName()}.${chromosome}"
        """
}

process indexGenotypes {
    label 'samtools'
    label 'mediumMemory'

    tag "${vcfFile}"

    input:
        path vcfFile
    output:
        tuple path("${vcfFile}"), path("${vcfFile}.tbi")
    script:
        """
        tabix -p vcf ${vcfFile}
        """
}

process concatenateAndSort {
    label 'bcftools'
    label 'mediumMemory'

    input:
        path vcfFiles
    output:
        path "${params.cohortName}-aligned.vcf.gz"
    script:
        """
        bcftools concat -Oz -o ${params.cohortName}-aligned.vcf.gz ${vcfFiles}
        """
}
