#!/usr/bin/env nextflow
/*
 *  REMOVE REALLY LOW QUALITY SAMPLES AND SNVS
 *  ==========================================
 *
 *  (...and align the genotypes to a human reference sequence)
 *
 *  This script performs basic quality control on our cohort data, 
 *  and should be done before sample quality control and snv quality
 *  control. We remove the samples and snvs that are obviously really
 *  bad and do not require careful treatment. 
 *
 *  We begin by selecting all variants in the input cohort data that
 *  are duplicated and then removing them from the data. Next we
 *  extract the cohort genotypes in the cohort data to a vcf file 
 *  and then align these genotpyes to the reference sequence. We then
 *  select only snvs that are biallelic with respect to this dataset 
 *  (i.e. exactly two alleles for a given snp can be found when 
 *  looking at all samples in the cohort). We then rebuild the cohort
 *  data (the plink bed, bim, and fam files) using the input fam file
 *  and this new aligned and filtered cohort genotypes vcf file. 
 *
 *  Finally, we perform a set of hard filtering transformations to 
 *  the rebuilt cohort data, which includes removing snvs with high 
 *  missingness rates, or those seriously far from hardy weinberg
 *  equilibrium and so on. 
 *
 *******************************************************************/

nextflow.enable.dsl=2

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/base.nf"

include {
    checkInputParams;
    getInputChannels;
    selectDuplicatedVariants;
    removeDuplicatedVariants;
    alignGenotypesToReference;
    selectBiallelicSnvs;
    rebuildCohortData;
    removeReallyLowQualitySamplesAndSnvs;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/basicQualityControl.nf"


workflow {

    checkInputParams()

    (cohortData,
     referenceSequence) = getInputChannels()

    duplicatedVariantIds \
        = selectDuplicatedVariants(
            cohortData)

    filteredCohortData \
        = removeDuplicatedVariants(
            cohortData,
            duplicatedVariantIds)

    alignedGenotypeSet \
        = alignGenotypesToReference(
            filteredCohortData,
            referenceSequence)

    biallelicGenotypeSet \
        = selectBiallelicSnvs(
            alignedGenotypeSet)

    alignedCohortData \
        = rebuildCohortData(
            biallelicGenotypeSet,
            filteredCohortData)

    basicFilteredCohortData \
        = removeReallyLowQualitySamplesAndSnvs(
            alignedCohortData)
}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}
