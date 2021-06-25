#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {

    plinkBinaryFileset = getPlinkBinaryFileset().toList()
    chromosome = getChromosomes()
    geneticMap = getGeneticMapFiles()
    thousandGenomesReference = getThousandGenomesReference()
    palindromes = getPalindromicSnvs( plinkBinaryFileset )

    perChromosomeVcfFiles = convertPlinkBinaryToVcf( plinkBinaryFileset, palindromes, chromosome )
    perChromosomeVcfFiles
	.map { chrom, vcfFile -> tuple( chrom, vcfFile ) }
	.join( thousandGenomesReference )
	.join( geneticMap )
	.set { checkstrand_inputs }

    alignedVcfs = alignGenotypesToReference( checkstrand_inputs )
    alignedVcfs
        .map { chrom, vcfFile -> tuple( chrom, vcfFile ) }
        .join( thousandGenomesReference )
        .join( geneticMap )
        .set { phasing_inputs }

    phasedVcfFiles = phaseGenotypes( phasing_inputs )

/*
    prePhasingQualityReports = getCheckStrandReports()
*/
	
}


def getPlinkBinaryFileset() {
    return channel.fromPath( params.outputDir + params.cohortName + '-clean.{bed,bim,fam}' )
}

def getChromosomes() {
    return channel.of(1..22, 'X')
}

def getThousandGenomesReference() {
    return channel.fromFilePairs( params.thousandGenomesDir + 'chr*1kg.phase3.v5a.vcf.gz', size: 1 )
	   	  .map { group_key, ref_file -> 
			 tuple( group_key.replaceFirst(/chr/,""), ref_file.first())
	    	  }
}

def getGeneticMapFiles() {
    return channel.fromFilePairs( params.thousandGenomesDir + 'plink.chr*.GRCh37.map', size: 1 )
		  .map { group_key, map_file ->
			 tuple( group_key.replaceFirst(/^plink\.chr/,""), map_file.first() )
		  }
}

process getPalindromicSnvs() {
    tag 'FILTER PALINDROMES'
    label 'Rbase'
    label 'mediumMemory'
    input:
        tuple path(bim), path(bed), path(fam)
    output:
        //publishDir path: "${params.outputDir}"
        path "${params.cohortName}.palindromic.snvs.txt"
    script:
        template "filterPalindromicVariants.r"
}

process convertPlinkBinaryToVcf() {
    tag 'BED+BIM+FAM ==> VCF'
    label 'plink2'
    label 'smallMemory'
    input:
	tuple path(bim), path(bed), path(fam)
	path palindromes
	val chromosome
    output:
	//publishDir path: "${params.outputDir}"
	tuple val("${chromosome}"), path("chr${chromosome}_forPhasing.vcf.gz")
    script:
	"""
	plink2 \
	    --bfile ${params.outputDir}${bed.baseName} \
	    --exclude ${palindromes} \
	    --export vcf-4.2 bgz id-paste='iid' \
	    --threads $task.cpus \
	    --chr ${chromosome} \
	    --out "chr${chromosome}_forPhasing"
	"""
}

process alignGenotypesToReference() {
    label 'beagle'
    input:
	tuple val(chromosome), path(vcfFile), path(refFile), path(geneticMap)
    output:
	//publishDir path: "${params.outputDir}"
        tuple val("${chromosome}"), path("chr${chromosome}-aligned.vcf.gz")
    script:
	"""
	conform-gt \
	    gt="${vcfFile}" \
	    ref="${refFile}" \
	    chrom=${chromosome} \
	    match=POS \
	    out="chr${chromosome}-aligned"
	"""
}

process phaseGenotypes() {
    label 'beagle'
    label 'bigMemory'
    input:
	tuple val(chromosome), path(vcfFile), path(refFile), path(geneticMap)
    output:
	publishDir path: "${params.outputDir}", mode: 'copy'
	path "chr${chromosome}-phased.vcf.gz"
    script:
	"""
	beagle \
	    gt="${vcfFile}" \
	    ref="${refFile}" \
	    map="${geneticMap}" \
	    chrom=${chromosome} \
	    burnin=3 \
	    iterations=12 \
	    ne=20000 \
	    impute=false \
	    nthreads=${task.cpus} \
	    out="chr${chromosome}-phased"
	"""
}

