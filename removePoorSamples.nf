/*
 *  Remove poor samples (sample QC)
 *
 ********************************************************************/
nextflow.enable.dsl=2

workflow {

    cohortData = getInputChannels()



}

def getInputChannels() {

    bed = channel.fromPath(
        "${params.inputDir}${params.cohortName}.bed")
    bim = channel.fromPath(
        "${params.inputDir}${params.cohortName}.bim")
    fam = channel.fromPath(
        "${params.inputDir}${params.cohortName}.fam")

    return bed.combine(bim).combine(fam)
}

process identifyIndivDiscSexinfo {

    memory plink_mem_req
    publishDir "${params.output_dir}/IndivDiscSexinfo", overwrite:true, mode:'copy'

    input:
        path bed
        path bim
        path fam

        path orig 
        path dup
        path lmiss
        path imiss

    output:
        path logfile
        path sexcheck_report
        path "${base}.hwe" 
        path imiss
        path lmiss

    script:
        base = bed.baseName
        logfile= "${base}.badsex"
        sexcheck_report = "${base}.sexcheck"
        imiss  = "${base}.imiss"
        lmiss  = "${base}.lmiss"
        """
        plink \
            --keep-allele-order \
            --bfile $base \
            --hardy \
            --check-sex $f_hi_female $f_lo_male \
            --missing \
            --out $base
        """
}