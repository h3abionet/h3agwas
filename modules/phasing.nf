include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/intensityPlot.nf"

include {
    checkCohortName;
    checkReferencePanelsDir;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
    getCohortData;
} from "${projectDir}/modules/base.nf"


def checkInputParams() {

    checkCohortName()
    checkReferencePanelsDir()
}

def getInputChannels() {
    return [
        getCohortData('clean'),
        getReferencePanels(),
        getGeneticMaps()]
}

process selectBiallelicSnvsWithBcftools {
    label 'bcftools'
    label 'mediumMemory'

    tag "referencePanel_chr${chromosome}"

    input:
    tuple val(chromosome), path(referencePanel), path(referencePanelIndex)
    output:
    tuple val(chromosome), path("chr${chromosome}.refpanel.biallelic.vcf.gz")
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

process selectGenotypeSetWithPlink {
    label 'plink2'
    label 'smallMemory'

    tag "cohortData"

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

process alignWithConformGt {
    label 'beagle'
    label 'bigMemory'

    tag "genotypeSet, referencePanel_chr${chromosome}"

    input:
        tuple path(cohortGenotypes), val(chromosome), path(referencePanel)
    output:
        tuple val(chromosome), path("${params.cohortName}.${chromosome}.vcf.gz")
    script:
        """
        conform-gt \
            ref=${referencePanel} \
            gt=${cohortGenotypes} \
            chrom=${chromosome} \
            out=${params.cohortName}.${chromosome}
        """
}

process phaseWithBeagle {
    label 'beagle'
    label 'bigMemory'

    tag "genotypeSubset_chr${chromosome}, referencePanel_chr${chromosome}, geneticMap_chr${chromosome}"

    input:
        tuple val(chromosome), path(alignedGenotypes), path(referencePanel), path(geneticMap)
    output:
        tuple val(chromosome), path("${params.cohortName}.${chromosomeString}.phased.vcf.gz")
    script:
        chromosomeString = ( chromosome < 10 ) ? "0${chromosome}" : "${chromosome}"
        """
        beagle \
            ref=${referencePanel} \
            map=${geneticMap} \
            gt=${alignedGenotypes} \
            chrom=${chromosome} \
            nthreads=${task.cpus} \
            window=${params.phasingWindowSize} \
            out=${params.cohortName}.${chromosomeString}.phased
        """
}

process indexWithTabix {
    label 'htslib'
    label 'mediumMemory'

    tag "haplotypeSubset_chr${chromosome}"

    input:
        tuple val(chromosome), path(vcfFile)
    output:
        tuple path("${vcfFile}"), path("${vcfFile}.tbi")
    script:
        """
        tabix -p vcf ${vcfFile}
        """
}

process concatenateWithBcftools {
    label 'bcftools'
    label 'bigMemory'

    tag "all haplotypeSubsets"

    input:
        path vcfFiles
    output:
        path "${params.cohortName}-phased.vcf.gz"
    script:
        """
        bcftools concat --threads $task.cpus -Oz -o ${params.cohortName}-phased.vcf.gz *.vcf.gz
        """
}

process rebuildCohortDataWithPlink() {
    label 'bigMemory'
    label 'plink2'

    tag "haplotypeSet, inputFam"

    input:
        path haplotypes
        tuple path(unphasedBed), path(unphasedBim), path(unphasedFam)
    output:
        publishDir path: "${params.outputDir}phasing", mode: 'copy'
        path "${params.cohortName}-phased.{bed,bim,fam,log}"
    script:
        """
        plink2 \
            --vcf ${haplotypes} \
            --fam ${unphasedFam} \
            --threads $task.cpus \
            --make-bed \
            --double-id \
            --out ${params.cohortName}-phased
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


def getReferencePanels() {
    return channel
        .of(1..22)
        .map{ it -> [
            it,
            file(params.referencePanelsDir + '*chr' + it + '.*.vcf.gz')[0],
            file(params.referencePanelsDir + '*chr' + it + '.*.vcf.gz.tbi')[0]]}
}

def getGeneticMaps() {
    return channel
        .of(1..22)
        .map{ it -> [
            it,
            file(params.geneticMapsDir + '*chr' + it + '.*.map')[0]]}
}
