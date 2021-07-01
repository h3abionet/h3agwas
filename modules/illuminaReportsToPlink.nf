include {
    checkCohortName;
    checkSampleReport;
    checkSnpReport;
    userEmailAddressIsProvided;
    checkEmailAdressProvided;
    getBasicEmailSubject;
    getBasicEmailMessage;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
        checkCohortName()
        checkSampleReport()
        checkSnpReport()
        checkEmailAdressProvided()
}

def getChunksFromGenotypeReports() {
    return channel
        .fromPath( params.inputDir + params.genotypeReport )
      	.splitText(
            by: 5000000,
            keepHeader: false,
            file: true,
            compress: false )
}

def getGenotypeReports() {

   return channel
            .fromPath( params.inputDir + "*_gtReport_*" )
}

def concatenateLgenChunks(lgenChunks) {
    return lgenChunks
        .collectFile(
            name: "${params.cohortName}.lgen",
            sort: true)
}

def getSampleReport() {
    return channel
                .fromPath( params.inputDir + params.sampleReport )
}

def getSnpReport() {
    return channel
                .fromPath( params.inputDir + params.snpReport )
}

process convertGenotypeReportsToLgen {
    tag "${chunks.baseName}"
    label 'smallMemory'
    label 'perl'
    cache 'lenient'
    input:
    	path chunks
    output:
    	path "*.lgen"
    script:
<<<<<<< HEAD
        template 'convertGenotypeReportsToLgen.py'
=======
        template 'convertGenotypeReportsToLgen.pl'
>>>>>>> 248854c9e6270e86ad730a982eaaadd40e21e82e

}

process concatenateLgenFiles() {

   input:
      path lgenFiles

   output:
      //publishDir path: "${params.outputDir}", mode: 'copy'
      path "${params.cohortName}.lgen"

   script:
      """
      cat ${lgenFiles} > "${params.cohortName}.lgen"
      """
}

process getMapFileFromSnpReport() {
    tag "$snpReport"
    label 'smallMemory'
    input:
        path snpReport
    output:
        path "${params.cohortName}.map"
    script:
        """
        cut \
            -f2-4 \
            -d',' ${snpReport} | \
        sed 's/,/ /g' | \
        awk '{print \$2,\$1,"0",\$3}' | \
        awk '\$1!="0"' | \
        sed '1d' > "${params.cohortName}.map"
        """
}

process getFamFileFromSampleReport() {
    tag "$sampleReport"
    label 'smallMemory'
    input:
        path sampleReport
    output:
        path "${params.cohortName}.fam"
    script:
        """
        grep \
             -v "Failed Sample" ${sampleReport} | \
        cut \
             -f2,13-14 \
            -d',' | \
        awk 'FS="," {print \$1,\$1,"0","0","-9","-9"}' | \
        sed '1d' > "${params.cohortName}.fam"
        """
}

process convertPlinkLongFormatToPlinkBinary() {
    tag "LGEN+MAP+FAM ==> BED+BIM+FAM"
    label 'mediumMemory'
    label 'plink'
    cache 'lenient'
    input:
        path "${params.cohortName}.lgen"
        path "${params.cohortName}.map"
        path "${params.cohortName}.fam"
    output:
        path "${params.cohortName}.{bed,bim,fam}"
    script:
        """
        plink \
            --lgen ${params.cohortName}.lgen \
            --map ${params.cohortName}.map \
            --fam ${params.cohortName}.fam \
	        --no-parents \
	        --no-sex \
	        --no-pheno \
	        --threads $task.cpus \
            --make-bed \
            --out ${params.cohortName}
        """
}

process alignGenotypeDataToReference() {
    label 'mediumMemory'
    label 'plink2'
    cache 'lenient'
    input:
	    path plinkBinaryFileset
	    path famFile
    output:
        path "temporary.vcf.gz"
    script:
    	plinkBase = famFile.baseName
        """
        plink2 \
            --bfile ${plinkBase} \
            --fa "${params.reference}" \
            --ref-from-fa force \
            --normalize \
            --threads $task.cpus \
            --export vcf-4.2 id-paste=iid bgz \
            --real-ref-alleles \
            --out temporary
        """
}

process filterSitesWithoutRefOrAltAlleles() {
    tag "BCFTOOLS VIEW"
    label 'mediumMemory'
    label 'bcftools'
    input:
        path tempVcfFile
    output:
        path "${params.cohortName}.vcf.gz"
    script:
        """
        bcftools \
            view \
            -m2 \
            --threads $task.cpus \
            -Oz \
            -o ${params.cohortName}.vcf.gz \
	    ${tempVcfFile}
        """
}

process getFinalPlinkBinaryFileset() {
    tag "VCF ==> PLINKBINARY"
    label 'mediumMemory'
    label 'plink2'
    input:
        path vcfFile
    output:
        publishDir path: "${params.outputDir}", mode: 'copy'
        path "${params.cohortName}-clean.{bed,bim,fam,log}"
    script:
        """
        plink2 \
            --vcf ${vcfFile} \
            --threads $task.cpus \
            --make-bed \
	    --double-id \
            --out ${params.cohortName}-clean
        """
}

def sendWorkflowExitEmail() {
    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: getBasicEmailSubject(),
<<<<<<< HEAD
            body: getBasicEmailMessage(),
            attach: [
                "${params.outputDir}manhattan.pdf",
                "${params.outputDir}qqplot.pdf"])
   }
}

=======
            body: getBasicEmailMessage())
   }
}
>>>>>>> 248854c9e6270e86ad730a982eaaadd40e21e82e
