/*
 *  PHASE COHORT GENOTYPE DATA
 *  ==========================
 *
 *
 ********************************************************************/
nextflow.enable.dsl=2

workflow {

/*
    test = channel.fromPath(
        params.inputDir + "test.29May21.d6d.vcf.gz")
    ref = channel.fromPath(
        params.inputDir + "ref.29May21.d6d.vcf.gz")
    target = channel.fromPath(
        params.inputDir + "target.29May21.d6d.vcf.gz")
*/

    (cohortData, legendFile) = getInputChannels()

    cohortVcfFile = convertPlinkBinaryToVcf(cohortData)
    cohortFrequencyFile = getFrequencyFile(cohortData)
    excludeFiles = checkStrand(
        cohortData,
        cohortFrequencyFile,
        legendFile)

    concatenatedExcludeFile = concatenate(excludeFiles) | view()
    excludeFile = removeDuplicateLines(concatenatedExcludeFile)

    cohortDataFiltered = removeUnalignedSnvs(
        cohortData,
        excludeFile) | view()

/*

    testBeagleWithGtInput(test)

    testBeagleWithGtAndRefInputs(target, ref)

    bref3 = makeBref3FromRef(ref)

    testBeagleWithGtAndBref3Inputs(target, bref3)

*/
}

def getInputChannels() {
    bed = channel.fromPath(
        params.inputDir + params.cohortName + ".bed")
    bim = channel.fromPath(
        params.inputDir + params.cohortName + ".bim")
    fam = channel.fromPath(
        params.inputDir + params.cohortName + ".fam")

    legendFile = channel.fromPath(
        params.databasesDir
        + 'Human/GRCh37/1000GP_Phase3/1000GP_Phase3_chr*.legend.gz')

    return [
        bed.combine(bim).combine(fam),
        legendFile]
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

process getFrequencyFile {
    label 'plink'
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        path "${cohortBed.getBaseName()}.frq"

    script:
        """
        plink \
            --bfile ${cohortBed.getBaseName()} \
            --freq \
            --out ${cohortBed.getBaseName()}
        """
}

process checkStrand {
    label 'checkStrand'
    label 'smallMemory'

    tag "${legendFile.getSimpleName()}"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
        path cohortFrequencyFile
        each path(legendFile)

    output:
        path "Exclude-${cohortBed.getBaseName()}-1000G.txt"

    script:
        """
        HRC-1000G-check-bim.pl \
            -b ${cohortBim} \
            -f ${cohortFrequencyFile} \
            -r ${legendFile} \
            -g \
            -p "AFR"
        """
}

def concatenate(inputFiles) {
    return inputFiles.collectFile(
        name: 'concatenated.txt',
        newLine: true)
}

process removeDuplicateLines {
    label 'smallMemory'

    input:
        path inputFile

    output:
        path "${outputFile}"

    script:
        outputFile = \
            "${inputFile.getBaseName()}" \
            + ".duplicatesRemoved." \
            + "${inputFile.getExtension()}"
        """
        sort ${inputFile} | sed '/^\$/d' | uniq > ${outputFile}
        echo '' >> ${outputFile}
        """
}

process removeUnalignedSnvs {
    label 'plink'
    label 'smallMemory'

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)
        path snvsToExclude

    output:
        path "${outputBase}.{bed, bim, fam}"
        
    script:
        outputBase = "${cohortBed.getBaseName()}.filtered"
        """
        plink \
            --bfile ${cohortBed.getBaseName()} \
            --exclude ${snvsToExclude} \
            --make-bed \
            --out ${outputBase}
        """
}

process testBeagleWithGtInput {
    label 'beagle'
    label 'smallMemory'
    publishDir "${params.outputDir}phased", mode: 'copy'

    input:
        path test

    output:
        path "out.gt.vcf.gz"

    script:
        """
        beagle gt=${test} out=out.gt
        """
}

process testBeagleWithGtAndRefInputs {
    label 'beagle'
    label 'smallMemory'
    publishDir "${params.outputDir}phased", mode: 'copy'

    input:
        path target
        path ref

    output:
        path "out.usingRef.gt.vcf.gz"

    script:
        """
        beagle \
            ref=${ref} \
            gt=${target} \
            out=out.usingRef.gt
        """

}

process makeBref3FromRef {
    label 'beagle'
    label 'smallMemory'

    input:
        path ref

    output:
        path "${output}"

    script:
        output = "${ref}".replaceFirst(/.vcf.gz/, ".bref3")
        """
        bref3 ${ref} > ${output}
        """
}

process testBeagleWithGtAndBref3Inputs {
    label 'beagle'
    label 'smallMemory'
    publishDir "${params.outputDir}/phased", mode: 'copy'

    input:
        path target
        path bref3

    output:
        path "out.bref3.vcf.gz"

    script:
        """
        beagle ref=${bref3} gt=${target} out=out.bref3
        """

}
