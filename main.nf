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
include { split_vcf} from "./modules/vcf.nf"
include { convertvcfin } from "./formatdata/format_vcf.nf"
include {qc_dup} from './qc/workflow.nf'
include {assoc} from './assoc/workflow.nf'
include {imputation} from './imputation/workflow.nf'
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
  data=null
  vcf_qc = null
  bimbam = null
  bgen = null
  bgenlist = null
  bgen_sample = null
  bfile=null
  build_cur=params.build_genome
  vcf=null
  type_impute=null
  if(params.qc_topbottom==1 || params.qc_topbottom){

  }
  if (params.qc_dup == 1 || params.qc_dup) {
        qc_dup(data, plink_qc, "${params.output_dir}/dup/")
        plink_qc=qc_dup.out.plink
        data=qc_dup.out.data
  }
  if (params.qc == 1 || params.qc) {
        qc(data, plink_qc, "${params.output_dir}/qc/")
        plink_qc=qc.out.plink
  }
  if (params.qc_michigan == 1 || params.qc_michigan) {
       qc_michigan(plink_qc)
       plink_qc = qc_michigan.out.plink
  }
  if (params.convertinvcf==1 || params.convertinvcf){
     format_plink_invcf(plink_qc, "${params.output_dir}/vcf/", params.output) 
     vcf_qc= format_plink_invcf.out.vcf
  }
  if(params.vcf_split_chro){
    split_vcf(vcf_qc,"${params.output_dir}/vcf/split", params.output)
    vcf_qc=split_vcf.out.vcf
  }
  if (params.vcf_convertbetwen_build==1 || params.vcf_convertbetwen_build){
    crossmap_vcf(vcf_qc) 
    build_cur=params.build_genome_convert
    vcf_qc = crossmap_vcf.out.vcf
 }
  if (params.impute==1 || params.impute){
     imputation(vcf_qc,build_cur, type_impute)
     vcf_qc = imputation.out.vcf
     type_impute= imputation.out.type_impute 
  }
 balise_convertvcf=(params.convertvcfinplink==1 || params.convertvcfinbimbam==1)
 if(balise_convertvcf){
   convertvcfin(vcf_qc, build_cur,"${params.output_dir}/convertvcf", params.output, type_impute)
   bfile = convertvcfin.out.plink
   vcf = convertvcfin.out.vcf
 }
 if(params.association==1 | params.assoc==1){
   assoc(data,bfile, vcf, bimbam, bgen,bgenlist, bgen_sample)
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
