
def getGenotypeReports() {
  
   return channel
            .fromPath( params.dataDir + "*_gtReport_*" )

}

def getSampleReport() {

   return channel
            .fromPath( params.dataDir + params.sampleReport )

}

def getSnpReport() {

   return channel
            .fromPath( params.dataDir + params.snpReport )

}



process convertGenotypeReportsToLgen() {

   label 'bigMemory'
   tag "${genotypeReports.baseName}"

   input:
      path genotypeReports

   output:
      publishDir path: "${params.outDir}", mode: 'copy'
      path "*.lgen"

   script:

      gsgt2lgenTemplate = "${projectDir}/modules/templates/gsgt2lgen.py"

      """
      python \
        ${gsgt2lgenTemplate} \
	     ${params.dataDir}${genotypeReports.baseName}.gz
      """

}

process concatenateLgenFiles() {

   input:
      path lgenFiles

   output:
      publishDir path: "${params.outDir}", mode: 'copy'
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
      publishDir path: "${params.outDir}", mode: 'copy'
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
      publishDir path: "${params.outDir}", mode: 'copy'
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

   input:
      path "${params.cohortName}_lgen"
      path "${params.cohortName}_map"
      path "${params.cohortName}_fam"

   output:
      publishDir path: "${params.outDir}", mode: 'copy'
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