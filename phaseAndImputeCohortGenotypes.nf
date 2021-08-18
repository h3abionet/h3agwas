/*
 *  PHASE COHORT GENOTYPE DATA
 *  ==========================
 *
 *  This script estimates the underlying haplotypes of the genotypes
 *  provided by the user, and imputes them, using reference panels and 
 *  genetic maps. The input is a plink binary trio (bed, bim, fam) 
 *  data set, below colectively referred to as `cohortData`. 
 *
 *  First the cohort data is read into a channel, along with the
 *  reference panels and genetic maps. We assume there is one panel
 *  and one map per chromosome. Then the reference panels are filtered
 *  such that only biallelic snvs are used. The genotype set of the
 *  cohort are extracted and saved temporarily in the VCF format.    
 *
 *  The genotype set of the cohort is then aligned to the reference 
 *  panels. Any snvs in this genotype set not found in the reference
 *  panels are removed. When aligning the genotype set to a given 
 *  reference panel, only the snvs on the chromosome represented by 
 *  the panel are output. Therefore a set of genotype subsets are 
 *  produced by the aligning step, where each subset contains snvs
 *  on a given chromosome.
 *
 *  Each of these genotype subsets is then phased with beagle, using
 *  the corresponding reference panel and genetic map. The snvs from 
 *  the panel are also imputed into the cohort data at this stage. 
 *  the result of this phasing step are imputed haplotype subsets. 
 *
 *  Finally, these haplotype subsets are indexed and then concatenated
 *  together to form a single haplotype set for the cohort, in vcf
 *  format. Using this file and the input cohort fam file we rebuild
 *  the plink binary trio (bed, bim, fam) for the phased cohort data.    
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    printWorkflowExitMessage;
} from "${projectDir}/modules/base.nf"

include {
    checkInputParams;
    getInputChannels;
    selectBiallelicSnvsWithBcftools;
    selectAutosomalGenotypeSet;
    recodeMaleXHaploidAsDiploid;
    alignWithConformGt;
    phaseWithBeagle;
    indexWithTabix;
    concatenateWithBcftools;
    rebuildCohortDataWithPlink;
    sendWorkflowExitEmail;
} from "${projectDir}/modules/phasing.nf"


workflow {

    checkInputParams()

    (inputCohortData,
     referencePanels,
     geneticMaps) = getInputChannels()

    filteredReferencePanels = selectBiallelicSnvsWithBcftools(
        referencePanels)

    genotypeSet = selectAutosomalGenotypeSet(
        inputCohortData)

    alignedGenotypeSubsets = alignWithConformGt(
        genotypeSet.combine(filteredReferencePanels))

    haplotypeSubsets = phaseWithBeagle(
        alignedGenotypeSubsets
            .join(filteredReferencePanels)
            .join(geneticMaps))

    indexedHaplotypeSubsets = indexWithTabix(
        haplotypeSubsets)

    haplotypeSet = concatenateWithBcftools(
        indexedHaplotypeSubsets.collect())

    phasedCohortData = rebuildCohortDataWithPlink(
        haplotypeSet, inputCohortData)
}

workflow.onComplete {
    printWorkflowExitMessage()
    sendWorkflowExitEmail()
}

