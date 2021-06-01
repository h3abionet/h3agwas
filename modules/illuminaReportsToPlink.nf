
def getGenotypeReports() {
  
   return channel
            .fromPath( params.inputDir + "*_gtReport_*" )

}

def getSampleReport() {

   return channel
            .fromPath( params.inputDir + params.sampleReport )

}

def getSnpReport() {

   return channel
            .fromPath( params.inputDir + params.snpReport )

}



process convertGenotypeReportsToLgen() {

   label 'bigMemory'
   label 'datatable'
   tag "${genotypeReports.baseName}"

   input:
      path genotypeReports

   output:
      publishDir path: "${params.outputDir}", mode: 'copy'
      path "*.lgen"

   script:

      gsgt2lgenTemplate = "${projectDir}/modules/templates/gsgt2lgen.py"

      """
      python \
        ${gsgt2lgenTemplate} \
	     ${params.inputDir}${genotypeReports.baseName}.gz \
        ${params.threads}
      """

}

process concatenateLgenFiles() {

   input:
      path lgenFiles

   output:
      publishDir path: "${params.outputDir}", mode: 'copy'
      path "${params.cohortName}.lgen"

   script:
      """
      cat ${lgenFiles} > "${params.cohortName}.lgen"
      """
   
}

process getMapFileFromSnpReport() {

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

   label 'plink'

   input:
      path "${params.cohortName}_lgen"
      path "${params.cohortName}_map"
      path "${params.cohortName}_fam"

   output:
      publishDir path: "${params.outputDir}", mode: 'copy'
      path "${params.cohortName}.{bed,bim,fam,log}"

   script:
      """
      plink \
         --lfile "${params.cohortName}_lgen" \
         --map "${params.cohortName}_map" \
         --fam "${params.cohortName}_fam" \
         --make-bed \
         --out ${params.cohortName}
      """
}