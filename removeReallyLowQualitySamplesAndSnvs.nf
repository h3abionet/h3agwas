#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    getInputChannels;
    removeReallyLowQualitySamplesAndSnvs;
} from "${projectDir}/modules/basicQualityControl.nf"

workflow {

    cohortData = getInputChannels()

    basicFilteredCohortData \
        = removeReallyLowQualitySamplesAndSnvs(
            cohortData)
}