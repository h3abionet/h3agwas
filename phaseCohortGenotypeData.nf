/*
 *  PHASE COHORT GENOTYPE DATA
 *  ==========================
 *
 *
 ********************************************************************/

nextflow.enable.dsl=2

params.inputDir = "${launchDir}/input"
params.outputDir = "${launchDir}/output"

workflow {

    input = channel.fromPath(
        params.inputDir + "test.29May21.d6d.vcf.gz")

    testBeagleContainer(input)

}

process testBeagleContainer {
    container 'sickleinafrica/beagle'
    publishDir "${params.outputDir}/phased", mode: 'copy'

    input:
        path test

    output:
        path "out.gt.vcf.gz"

    script:
        """
        beagle gt=${test} out=out.gt
        """
}