include {extractPheno} from '../modules/phenotype.nf'
include {wf_prepare_pheno} from '../modules/phenotype.nf'
include {clean_plink;indep_pairwise} from '../modules/utils_plink.nf'
include {buildindex as buildindex_vcf} from '../modules/vcf.nf'
include {list_chro} from '../modules/utils_plink.nf'
include {compute_pcs_small} from '../modules/pcs.nf'
include {saige} from './process/saige.nf'
include {regenie} from './process/regenie.nf'
include {merge_sumstat} from './process/all.nf'
workflow check_params{
  take : 
    data
    bfile
    vcf
    bimbam
    bgen
    bgenlist
    bgen_sample
 main :
 if(data==null){
   data_ch=channel.fromPath(params.data, checkIfExists:true)
 }else{
  data_ch = data
 }
 balise_vcf=true
if (bfile==null){
  if(params.bfile!=''){
        inpat=params.bfile
  }else{
        inpat = "${params.input_dir}/${params.input_pat}"
   }
  bfile=Channel.fromPath("${inpat}.bed",checkIfExists:true).combine(Channel.fromPath("${inpat}.bim",checkIfExists:true)).combine(Channel.fromPath("${inpat}.fam",checkIfExists:true))

 }
 if(vcf==null){
   balise_vcf=false
   if(params.vcf_list!='') {
   balise_vcf=true
      vcf=Channel.fromPath(file(params.vcf_list).readLines(), checkIfExists:true)
    }    else if(params.vcf!=''){
     balise_vcf=true
     vcf=channel.fromPath(params.vcf, checkIfExists:true) 
   }
 }
 if(balise_vcf){
   vcf=buildindex_vcf(vcf)
 }
 balise_bgen=true
 balise_bgenlist=true
 if(bgenlist==null){
   balise_bgenlist=false
   if(params.bgen_list!='') {
      balise_bgenlist=true
      bgenlist=Channel.fromPath(file(params.bgen_list).readLines(), checkIfExists:true)
      bgen_sample=channel.fromPath(params.bgen_sample, checkIfExists:true)
    }
 }
if(bgen==null){
      balise_bgen=false
   if(params.bgen!='') {
      balise_bgen=true
      bgen=Channel.fromPath(params.bgen, checkIfExists:true)
      bgen_sample=channel.fromPath(params.bgen_sample, checkIfExists:true)
    }
 }
 extractPheno(data_ch, channel.from(params.pheno), channel.from(params.covariates), channel.from(params.gxe))
 if(params.snps_include_rel)incl_rel=channel.fromPath(params.snps_include_rel, checkIfExists:true) else incl_rel = channel.fromPath("01")
 if(params.snps_exclude_rel)excl_rel=channel.fromPath(params.snps_exclude_rel, checkIfExists:true) else excl_rel = channel.fromPath("02")
 clean_plink(bfile,incl_rel, excl_rel, extractPheno.out, channel.from(params.cut_maf_rel), channel.from(params.snps_thincount_rel), channel.from(params.snps_justacgt_rel))
 indep_pairwise(clean_plink.out, channel.of(params.plink_indep_pairwise), channel.of("${params.output_dir}/assoc/rel/"))
 
 /*computed pcs*/
compute_pcs_small(indep_pairwise.out, channel.from(params.add_pcs), channel.from("${params.output_dir}/assoc/pcs/"))
 /*plot pcs ?*/
 
 wf_prepare_pheno(extractPheno.out, channel.from(params.pheno), channel.from(params.covariates), channel.from(params.gxe),clean_plink.out, compute_pcs_small.out.eigen)
 pheno_type=params.phenotypes_type
 listchro=list_chro(bfile)
 listchro=listchro.flatMap{ list_str -> list_str.split() }
 bfile_rel=bfile
 if(params.sample_snps_rel==1){
  println('relatdness using snps independant')
  bfile_rel=indep_pairwise.out
 }else {
  println('full SNPs used for relatdness')
 }

 emit :
  data = wf_prepare_pheno.out.data
  pheno=  wf_prepare_pheno.out.pheno
  pheno_bin=  channel.from(params.phenotypes_type)
  covar = wf_prepare_pheno.out.covar
  covariates_type =channel.from(params.covariates_type)
  gxe = wf_prepare_pheno.out.gxe
  bfile = bfile
  bfile_rel = bfile_rel
  vcf_balise=balise_vcf
  vcf=vcf
  bgen=bgen
  bgenlist=bgenlist
  bgen_sample=bgen_sample
  bgen_balise=balise_bgen
  bgenlist_balise=balise_bgenlist
  listchro=listchro
}


workflow assoc {
 take :
   data
   bfile
   vcf
   bimbam
   bgen
   bgenlist
   bgen_sample
 main :
  println "In ASSOC top"
  check_params(data,bfile, vcf, bimbam, bgen,bgenlist,bgen_sample)
  if(params.saige){
   saige(check_params.out.data,   check_params.out.pheno,check_params.out.pheno_bin,check_params.out.covar, check_params.out.covariates_type, check_params.out.bfile,check_params.out.bfile_rel,check_params.out.vcf, check_params.out.vcf_balise,check_params.out.bgen, check_params.out.bgenlist,check_params.out.bgen_sample,check_params.out.bgen_balise, check_params.out.bgenlist_balise, check_params.out.listchro)
  }
if(params.regenie){
   regenie(check_params.out.data,   check_params.out.pheno,check_params.out.pheno_bin,check_params.out.covar, check_params.out.covariates_type, check_params.out.bfile,check_params.out.bfile_rel,check_params.out.vcf, check_params.out.vcf_balise,check_params.out.bgen, check_params.out.bgenlist,check_params.out.bgen_sample,check_params.out.bgen_balise, check_params.out.bgenlist_balise, check_params.out.listchro)
 }
 // add pcs
       println "In ASSOC middle"
}
