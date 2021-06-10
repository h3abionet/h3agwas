include {
    checkCohortName;
    checkSampleReport;
    checkSnpReport;
} from "${projectDir}/modules/base.nf"

def checkInputParams() {
    checkCohortName()
    checkSampleReport()
    checkSnpReport()
}

def getChunksFromGenotypeReports() { 
   return channel
      	     .fromPath( params.inputDir + "*_gtReport_*" )
      	     .splitText( by: 5000000,
             	         keepHeader: false,
                  	 file: true,
                  	 compress: false )
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
    cache 'lenient'
    input:
    	path chunks
    output:
    	publishDir path: "${params.outputDir}"
    	path "*.lgen"
    script:
    	template "convertGenotypeReportsToLgen.pl"
}

process getMapFileFromSnpReport() {
   tag "$snpReport"
   label 'smallMemory'
   input:
      path snpReport
   output:
      publishDir path: "${params.outputDir}", mode: 'copy'
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
      publishDir path: "${params.outputDir}", mode: 'copy'
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
   input:
      path "${params.cohortName}.lgen"
      path "${params.cohortName}.map"
      path "${params.cohortName}.fam"
   output:
      publishDir path: "${params.outputDir}", mode: 'copy'
      path "${params.cohortName}.{bed,bim,fam,log}"
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

