/*
 *  PHASE COHORT GENOTYPE DATA
 *  ==========================
 *
 *
 ********************************************************************/
nextflow.enable.dsl=2

workflow {

    test = channel.fromPath(
        params.inputDir + "test.29May21.d6d.vcf.gz")
    ref = channel.fromPath(
        params.inputDir + "ref.29May21.d6d.vcf.gz")
    target = channel.fromPath(
        params.inputDir + "target.29May21.d6d.vcf.gz")

    cohortData = getInputChannels() | view()

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

    return bed.combine(bim).combine(fam)
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
