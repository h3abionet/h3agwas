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

def getInputChannels() {
    return [
        getGenotypeReports(),
        getSampleReport(),
        getSnpReport()]
}

def getGenotypeReports() { 
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


process convertGenotypeReportsToLgen() {
    label 'smallMemory'
    //label 'datatable'
    tag "${genotypeReports.baseName}"
    input:
        path genotypeReportChunk
    output:
        publishDir path: "${params.outputDir}"
        path "*.lgen"
    script:
        //template 'convertGenotypeReportsToLgen.pl'
        """
        perl \
            "${launchDir}/templates/convertGenotypeReportsToLgen.pl" \
            ${genotypeReportChunk} \
            ${params.numberOfGtReportHeaderLines}
        """
}

/*
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
*/

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
      path "${params.cohortName}.lgen"
      path "${params.cohortName}.map"
      path "${params.cohortName}.fam"
   output:
      publishDir path: "${params.outputDir}", mode: 'copy'
      path "${params.cohortName}.{bed,bim,fam,log}"
   script:
      """
      plink \
         --lfile "${params.cohortName}.lgen" \
         --map "${params.cohortName}.map" \
         --fam "${params.cohortName}.fam" \
         --make-bed \
         --out ${params.cohortName}
      """
}

def sendWorkflowExitEmail() {
    if (userEmailAddressIsProvided()) {
        sendMail(
            to: "${params.email}",
            subject: getBasicEmailSubject(),
            body: getBasicEmailMessage(),
            attach: [
                "${params.reportsDir}/report.html"])
   }
}