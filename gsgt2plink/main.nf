#!/usr/bin/env nextflow

nextflow.enable.dsl = 2






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
      path "*.lgen", emit: lgen_ch

   script:
      """
      python \
	     ${launchDir}/gsgt2lgen.py \
	     ${params.dataDir}${genotypeReports.baseName}.gz
      """

}

process concatenateLgenFiles() {

   input:
      path lgenFiles

   output:
      publishDir path: "${params.outDir}", mode: 'copy'
      path "gwas.lgen"

   script:
      """
      cat ${lgenFiles} > gwas.lgen
      """
   
}

process getMapFileFromSnpReport() {

   label 'smallMemory'

   input:
      path snpReport

   output:
      publishDir path: "${params.outDir}", mode: 'copy'
      path "gwas.map", emit: mapFile_ch

   script:
      """
      cut \
         -f2-4 \
         -d',' ${snpReport} | \
      sed 's/,/ /g' | \
      awk '{print \$2,\$1,"0",\$3}' | \
      awk '\$1!="0"' | \
      sed '1d' > gwas.map
      """

}

process getFamFileFromSampleReport() {

   label 'smallMemory'

   input:
      path sampleReport

   output:
      publishDir path: "${params.outDir}", mode: 'copy'
      path "gwas.fam", emit: famFile_ch

   script:
      """
      grep \
         -v "Failed Sample" ${sampleReport} | \
      cut \
         -f2,13-14 \
         -d',' | \
      awk 'FS="," {print \$1,\$1,"0","0","-9","-9"}' | \
      sed '1d' > gwas.fam
      """
}

process convertPlinkLongFormatToPlinkBinary() {

   input:
      path gwas_lgen
      path gwas_map
      path gwas_fam

   output:
      publishDir path: "${params.outDir}", mode: 'copy'
      path "cameroonScdGwasCohort.{bed,bim,fam,log}"

   script:
      """
      plink \
         --lfile ${gwas_lgen} \
         --map ${gwas_map} \
         --fam ${gwas_fam} \
         --make-bed \
         --out cameroonScdGwasCohort
      """
}


workflow {

   genotypeReports_ch = getGenotypeReports()
   sampleReport_ch = getSampleReport()
   snpReport_ch = getSnpReport()

   convertGenotypeReportsToLgen( genotypeReports_ch )
   getFamFileFromSampleReport( sampleReport_ch )
   getMapFileFromSnpReport( snpReport_ch )
   concatenateLgenFiles( convertGenotypeReportsToLgen.out.lgen_ch )
   convertPlinkLongFormatToPlinkBinary( convertGenotypeReportsToLgen.out.lgen_ch, 
                                        getFamFileFromSampleReport.out.famFile_ch, 
                                        getMapFileFromSnpReport.out.mapFile_ch )

}
