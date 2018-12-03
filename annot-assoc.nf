#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Shaun Aron
 *   	Rob Clucas
 *      Eugene de Beste
 *      Lerato Magosi
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description : pipeline annotation 
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths




def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin"]
// define param for 
annotation_model=["gemma","boltlmm", "plink", "head", "linear", "logistic", "fisher", "fastlmm", ""]
allowed_params+=annotation_model
annotation_param=["list_rs", "file_annotation", "file_gwas", "around_rs"]
allowed_params+=annotation_param
allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2"]
allowed_params+=allowed_params_head


params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}



def params_help = new LinkedHashMap(helps)


params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output_testing = "cleaned"
params.covariates = ""
outfname = params.output_testing


/* Defines the path where any scripts to be executed can be found.
 */

/**/
params.head_pval = "P_BOLT_LMM"
params.head_freq = "A1FREQ"
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.head_beta="BETA"
params.head_se="SE"
params.head_A1="ALLELE1"
params.head_A2="ALLELE0"

params.around_rs=100000
params.cut_maf = 0.01

params.loczm_bin  = ""
params.loczm_pop = "AFR"
params.loczm_build = "hg19"
params.loczm_source ="1000G_March2012"
params.loczm_gwascat = ""




/* Do permutation testing -- 0 for none, otherwise give number */
/*JT Append initialisation variable*/
params.bolt_impute2fidiid=""
/*gxe param : contains column of gxe*/
params.gxe=""


params.input_pat  = 'raw-GWA-data'

params.sexinfo_available = "false"


params.plink_mem_req = '750MB' // how much plink needs for this
params.other_process_memory = '750MB' // how much other processed need


plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_memory
max_plink_cores = params.max_plink_cores 
params.help = false


//data_ch = Channel.fromPath(params.data)

//---- Modification of variables for pipeline -------------------------------//


def getConfig = {
  all_files = workflow.configFiles.unique()
  text = ""
  all_files.each { fname ->
      base = fname.baseName
      curr = "\n\n*-subsection{*-protect*-url{$base}}@.@@.@*-footnotesize@.@*-begin{verbatim}"
      file(fname).eachLine { String line ->
	if (line.contains("secretKey")) { line = "secretKey='*******'" }
        if (line.contains("accessKey")) { line = "accessKey='*******'" }
        curr = curr + "@.@"+line 
      }
      curr = curr +"@.@*-end{verbatim}\n"
      text = text+curr
  }
  return text
}



// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}
if(params.boltlmm==1){
head_pval = "P_BOLT_LMM"
head_freq = "A1FREQ"
head_bp = "BP"
head_chr = "CHR"
head_rs = "SNP"
head_beta="BETA"
head_se="SE"
head_A1="ALLELE1"
head_A2="ALLELE0"
}else if(params.gemma==1){

}else{
head_pval=params.head_pval 
head_freq=params.head_freq
head_bp=params.head_bp
head_chr=params.head_chr
head_rs=params.head_rs
head_beta=params.head_beta
head_se=params.head_se
head_A1=params.head_A1
head_A2=params.head_A2
}





println "\nTesting rs : ${params.list_rs}\n"
println "\nTesting data : ${params.data}\n"
println "Testing for phenotypes  : ${params.pheno}\n"
println "Testing for gwas file : ${params.file_gwas}\n"
println "Using covariates        : ${params.covariates}\n"
println "Around rs : ${params.around_rs}\n\n"

rs_label_ch = Channel.from(params.list_rs.split(","))
gwas_ch = Channel.fromPath(params.file_gwas)

process ExtractInfoRs{
    input:
       file(gwas_file) from gwas_ch
    each rs from rs_label_ch
    output :
      set rs, file(out_locus_rs) into locuszoom_ch
      set rs, file(out_info_rs) into infors_rs
      set rs, file(out_gwas_rs) into infors_gwas
    script:
      out="sub"
      out_locus_rs=out+"_"+rs+"_around.stat"
      out_gwas_rs=out+"_"+rs+"_gwas.stat"
      out_info_rs=out+"_"+rs+"_info.stat"
      """
     an_extract_rs.py --inp_resgwas  $gwas_file --chro_header $head_chr --pos_header $head_bp --rs_header $head_rs --pval_header $head_pval --beta_header ${head_beta} --freq_header  $head_freq --maf ${params.cut_maf} --list_rs $rs --around_rs ${params.around_rs} --out_head $out
      """
}

if(params.loczm_gwascat!=""){
loczm_gwascat=" --gwas-cat ${params.loczm_gwascat}"
}else{
loczm_gwascat=""
}

process PlotLocusZoom{
    input : 
       set rs, file(filegwas) from locuszoom_ch
    publishDir "${params.output_dir}/$rs", overwrite:true, mode:'copy'
    output :
       file("out_$rs/*.svg") into report_rs_ch
    script :
       """
       ${params.loczm_bin} --epacts  $filegwas --delim tab --refsnp  $rs --flank ${params.around_rs} --pop ${params.loczm_pop} --build ${params.loczm_build} --source ${params.loczm_source} --gwas-cat whole-cat_significant-only --svg  -p out --no-date 
       """
}

fileannot_ch = Channel.fromPath(params.list_file_annot).mix(infors_rs)

process ExtractAnnotation{
      input :
        set file(infoannot),rs, file(file_rs) from fileannot_ch
      script :    
         """  
         """
}

process PlotByGenotype{


}

