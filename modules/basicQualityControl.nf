include {
    checkOutputDir;
    checkCohortName;
    checkEmailAdressProvided;
    userEmailAddressIsProvided;
    getCohortData;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

def getInputChannels() {
    return getCohortData('input')
}

process removeReallyLowQualitySamplesAndSnvs {
    label 'smallMemory' 

    tag "Removing QC Phase 1"

    input:
        tuple path(cohortBed), path(cohortBim), path(cohortFam)

    output:
        publishDir "${params.outputDir}", mode: 'copy'
        path "${params.cohortName}.basicFiltered.{bed,bim,fam}"

    script:
        """
        plink \
            --keep-allele-order \
            --autosome \
            --bfile ${cohortBed.getBaseName()} \
            --mind 0.1 \
            --geno 0.1 \
            --make-bed \
            --out temp1
        plink \
            --keep-allele-order \
            --bfile temp1  \
            --mind $params.cut_mind \
            --make-bed \
            --out temp2
        plink \
            --keep-allele-order \
            --bfile temp2 \
            --geno $params.cut_geno \
            --make-bed \
            --out temp3
        plink \
            --keep-allele-order \
            --bfile temp3 \
            --maf $params.cut_maf \
            --make-bed \
            --out temp4
        plink \
            --keep-allele-order \
            --bfile temp4  \
            --hwe $params.cut_hwe \
            --make-bed \
            --out ${params.cohortName}.basicFiltered     
        """
}