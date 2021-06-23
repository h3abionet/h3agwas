/*
 *  PHASE COHORT GENOTYPE DATA
 *  ==========================
 *
 *
 ********************************************************************/
nextflow.enable.dsl=2

workflow {

    cohortData = getInputChannels() | view()

    cohortGenotypes = convertPlinkBinaryToVcf(cohortData) | view()

    alignGenotypesToReference(cohortGenotypes) | view()
}

def getInputChannels() {
    bed = channel.fromPath(
        params.outputDir + params.cohortName + ".clean.bed")
    bim = channel.fromPath(
        params.outputDir + params.cohortName + ".clean.bim")
    fam = channel.fromPath(
        params.outputDir + params.cohortName + ".clean.fam")

    return bed.combine(bim).combine(fam)
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
            --export vcf \
            --out ${cohortBed.getBaseName()}
        gzip ${cohortBed.getBaseName()}.vcf
        """
}

process alignGenotypesToReference {
    label 'beagle'
    label 'smallMemory'

    input:
        path cohortGenotypes
    output:
        stdout
    script:
        """
        conform-gt        
        """
}
