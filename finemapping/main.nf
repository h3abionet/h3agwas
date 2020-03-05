#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2019
 *
 *
 * Description : pipeline to do a finemapping 
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}



def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts", "max_forks", "cut_maf", "phenotype", "accessKey", "access-key", "secretKey", "secret-key",  "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin"]
params_mf=[""]
allowed_params+=params_mf



def params_help = new LinkedHashMap(helps)

params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"

params.gcta_bin="gcta64"

## paramater
params.n_pop=10000

## params file input
params.head_pval = "P_BOLT_LMM"
params.head_freq = ""
params.head_n = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.head_beta="BETA"
params.head_se="SE"
params.head_A1="ALLELE0"
params.head_A2="ALLELE1"

params.cut_maf=0.01
params.plink_mem_req="6GB"

## gcta parameters
params.gcta_mem_req="15GB"
params.cojo_p=1e-7
params.cojo_wind=10000
params.gcta_cpus_req = 1
params.cojo_slct=1
params.cojo_slct_other=""
params.cojo_actual_geno=0
params.big_time='100h'

params.finemap_bin="finemap"
params.caviarbf_bin="caviarbf"
paintor.paintor_bin="PAINTOR"

params.chro=""
params.begin_seq=""
params.end_seq=""
if(params.chro=="" | params.begin_seq=="" | params.end_seq==""){
error('chro, begin_seq or end_seq not initialise')
}



params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}


bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()

raw_src_ch= Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


gwas_extract_plk=Channel.create()
plink_subplk=Channel.create()
raw_src_ch.separate( gwas_extract_plk, plink_cojo, plink_other) { a -> [ a, a, a] }


gwas_file=Channel.fromPath(params.file_gwas)
// plink 

process ExtractPosition{
  input :
     file(file) from gwas_file
     file(bed),file(bim),file(fam) from plk
  output :
    file("${gwas}.gcta") into gcta_gwas
    file("${gwas}_finemap.z") into  finemap_gwas
    file("${gwas}_caviar.z") into range_gwas
    file("${gwas}.range") into range_plink
  script :
    out=chro+"_"+params.begin_seq+" "+params.end_seq
    """
    fine_extract_sig.py --inp_resgwas $GWASTmp --chro ${params.chro} --begin ${params.begin_seq}  --end ${params.end_seq} --chro_header ${params.head_chr} --pos_header $params.head_bp} --beta_header ${parans.head_beta} --se_header ${params.head_se} --a1_header ${params.head_A1} --a2_header ${params.head_A2} --freq_header  ${params.head_freq} --bim_file  $bim --rs_header ${parans.head_rs} --out_head $out --p_header ${params.head_pval}  --n ${params.n_pop}
    """
}

process SubPlink{
 

}
