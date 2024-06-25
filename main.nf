#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *      Jean-Tristan Brandenburg
 *      Scott Hazelhurst
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2024
 *
 *
 * Description  : Nextflow pipeline for association .
 *
 *(C) University of the Witwatersrand, Johannesburg, 2016-2024
 *    on behalf of the H3ABioNet Consortium
 *This is licensed under the MIT Licence. See the "LICENSE" file for details
 */


// quality control
include { qc } from "./qc/workflow.nf"
include { qc_michigan } from "./qc/workflow.nf"
include { format_plink_invcf } from "./formatdata/workflow.nf"
include { checkresume } from "./modules/fct_groovy.nf"
include { crossmap_vcf } from "./convertdatabuild/workflow.nf"
include { mtag } from "./mtag/workflow.nf"
include { ldsc } from "./heritability/workflow.nf"
include { convertvcfin } from "./formatdata/format_vcf.nf"
//include { assoc} from "./assoc/assoc.nf"
//include { assoc} from "./assoc/assoc.nf"
//nextflow.enable.moduleBinaries = true

workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}

workflow {
  checkresume()
  //if( !nextflow.version.matches('21.04+') ) {
  //  println "This workflow requires Nextflow version 21.04 or greater -- You are running version $nextflow.version"
  //  exit 1
  //}
  plink_qc=null
  vcf_qc = null
  build_cur=params.build_genome
  if (params.qc == 1 || params.qc) {
        qc()
        plink_qc=qc.out.plink
  }
  if (params.qc_michigan == 1 || params.qc_michigan) {
       qc_michigan(plink_qc)
       plink_qc = qc_michigan.out.plink
  }
  if (params.convertinvcf==1 || params.convertinvcf){
     format_plink_invcf(plink_qc, "${params.output_dir}/plk_vcf", params.output) 
     vcf_qc= format_plink_invcf.out.vcf
  }
  if (params.vcf_convertbetwen_build==1 || params.vcf_convertbetwen_build){
   crossmap_vcf(vcf_qc) 
  build_cur=params.build_genome_convert
   vcf_qc = crossmap_vcf.out.vcf
  }
 balise_convertvcf=(params.convertvcfinplink==1 || params.convertvcfinbimbam==1)
 if(balise_convertvcf){
   convertvcfin(vcf_qc, build_cur,"${params.output_dir}/convertvcf", params.output)
 }
  sumstat=null
  if (params.mtag==1){
    mtag("$params.output_dir/mtag")
    sumstat = mtag.out.sumstat_clean
  }
  if(params.ldsc == 1 ){
   ldsc(sumstat, "$params.output_dir/heritability/ldsc")
  }
}
