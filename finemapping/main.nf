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

allowed_params = ["input_dir","input_pat","output","output_dir","data","covariates", "work_dir", "scripts", "max_forks", "cut_maf", "phenotype", "accessKey", "access-key", "secretKey", "secret-key",  "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin", "batch", "batch_col" ,"samplesize", "manifest", "region", "AMI", "queue", "strandreport"]
params_bin=["finemap_bin", "finemap_bin", "paintor_bin","plink_bin", "caviarbf_bin", "gcta_bin"]
params_mf=["chro", "begin_seq", "end_seq", "n_pop","threshold_p", "n_causal_snp"]
params_cojo=["cojo_slct_other", "cojo_top_snps","cojo_slct", "cojo_actual_geno"]
params_filegwas=[ "file_gwas", "head_beta", "head_se", "head_A1", "head_A2", "head_freq", "head_chr", "head_bp", "head_rs", "head_pval"]
params_memcpu=["gcta_mem_req","plink_mem_req", "other_mem_req","gcta_cpus_req", "fm_cpus_req"]
allowed_params+=params_mf
allowed_params+=params_cojo
allowed_params+=params_filegwas
allowed_params+=params_bin
allowed_params+=params_memcpu



def params_help = new LinkedHashMap(helps)

params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"

params.gcta_bin="gcta64"

// paramater
params.n_pop=10000

// params file input
params.head_pval = "P_BOLT_LMM"
params.head_freq = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.head_beta="BETA"
params.head_se="SE"
params.head_A1="ALLELE0"
params.head_A2="ALLELE1"

params.cut_maf=0.01
params.plink_mem_req="6GB"

// gcta parameters
params.gcta_mem_req="15GB"
params.cojo_p=1e-7
params.cojo_wind=10000
params.gcta_cpus_req = 1
params.fm_cpus_req = 5
params.fm_mem_req = "20G"
params.cojo_slct=1
params.cojo_slct_other=""
params.cojo_actual_geno=0
params.big_time='100h'

params.threshold_p=5*10**-8
params.n_causal_snp=3
params.caviarbf_avalue="0.1,0.2,0.4"
params.paintor_fileannot=""
params.paintor_annot=""


params.finemap_bin="finemap"
params.caviarbf_bin="caviarbf"
params.paintor_bin="PAINTOR"
params.plink_bin="plink"

params.chro=""
params.begin_seq=""
params.end_seq=""
if(params.file_gwas==""){
error('file_gwas option not initialise')
}
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
raw_src_ch.separate( gwas_extract_plk, plink_subplk) { a -> [ a, a] }


gwas_file=Channel.fromPath(params.file_gwas)
// plink 

process ExtractPositionGwas{
  input :
     file(filegwas) from gwas_file
     set file(bed),file(bim),file(fam) from gwas_extract_plk
  output :
    file("${out}.gcta") into gcta_gwas
    file("${out}_finemap.z") into  (finemap_gwas_cond, finemap_gwas_sss)
    file("${out}_caviar.z") into caviarbf_gwas
    file("${out}.paintor") into paintor_gwas
    file("${out}.range") into range_plink
  script :
    out=params.chro+"_"+params.begin_seq+"_"+params.end_seq
    """
    fine_extract_sig.py --inp_resgwas $filegwas --chro ${params.chro} --begin ${params.begin_seq}  --end ${params.end_seq} --chro_header ${params.head_chr} --pos_header ${params.head_bp} --beta_header ${params.head_beta} --se_header ${params.head_se} --a1_header ${params.head_A1} --a2_header ${params.head_A2} --freq_header  ${params.head_freq} --bim_file  $bim --rs_header ${params.head_rs} --out_head $out --p_header ${params.head_pval}  --n ${params.n_pop}
    """
}

process SubPlink{
  input :
     set file(bed),file(bim),file(fam) from plink_subplk
     file(range) from range_plink
  output :
     set file("${out}.bed"),file("${out}.bim"),file("${out}.fam") into (subplink_ld, subplink_gcta)
  script : 
     plk=bed.baseName
     out=plk+'_sub'
     """
     ${params.plink_bin} -bfile $plk  --keep-allele-order --extract  range  $range --make-bed -out  $out
     """
}

process ComputedLd{
   input : 
      set file(bed),file(bim),file(fam) from subplink_ld
   output :
       file("$outld") into (ld_fmcond, ld_fmsss,ld_caviarbf, ld_paintor)
   script :
    outld=params.chro+"_"+params.begin_seq+"_"+params.end_seq+".ld"
    plk=bed.baseName
    """
     ${params.plink_bin} --r2 square0 yes-really -bfile $plk -out "tmp"
    sed 's/\\t/ /g' tmp.ld > $outld
    """
}

process ComputedFineMapCond{
  cpus params.fm_cpus_req
  memory params.fm_mem_req
  input :
    file(ld) from ld_fmcond 
    file(filez) from finemap_gwas_cond
  publishDir "${params.output_dir}/fm_cond", overwrite:true, mode:'copy'
  output :
    file("${out}.snp") into res_dmcond
    set file("${out}.config"), file("${out}.cred"), file("${out}.log_cond")
  script:
  fileconfig="config"
  out=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_cond" 
  """ 
  echo "z;ld;snp;config;cred;log;n_samples" > $fileconfig
  echo "$filez;$ld;${out}.snp;${out}.config;${out}.cred;${out}.log;${params.n_pop}" >> $fileconfig
  ${params.finemap_bin} --cond --in-files $fileconfig   --log --cond-pvalue ${params.threshold_p}  --n-causal-snps ${params.n_causal_snp}
  """
}

process ComputedFineMapSSS{
  memory params.fm_mem_req
  cpus params.fm_cpus_req
  input :
    file(ld) from ld_fmsss
    file(filez) from finemap_gwas_sss
  publishDir "${params.output_dir}/fm_sss", overwrite:true, mode:'copy'
  output :
    file("${out}.snp") into res_fmsss
    set file("${out}.config"), file("${out}.cred${params.n_causal_snp}"), file("${out}.log_sss")
  script:
  fileconfig="config"
  out=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_sss"
  """
  echo "z;ld;snp;config;cred;log;n_samples" > $fileconfig
  echo "$filez;$ld;${out}.snp;${out}.config;${out}.cred;${out}.log;${params.n_pop}" >> $fileconfig
  ${params.finemap_bin} --sss --in-files $fileconfig  --n-threads ${params.fm_cpus_req}  --log --n-causal-snps ${params.n_causal_snp}
  """
}

process ComputedCaviarBF{
  memory params.fm_mem_req
  input :
    file(filez) from caviarbf_gwas
    file(ld) from ld_caviarbf
  publishDir "${params.output_dir}/caviarbf", overwrite:true, mode:'copy'
  output :
   file("$output") into res_caviarbf
  script :
   output=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_caviarbf"
   """
   ${params.caviarbf_bin} -z ${filez} -r $ld  -t 0 -a ${params.caviarbf_avalue} -c ${params.n_causal_snp} -o ${output} -n ${params.n_pop}
   """
}
if(params.paintor_fileannot=="" || params.paintor_annot==""){
paintor_fileannot=file('NOFILE')
}else{
paintor_fileannot=Channel.fromPath(params.paintor_fileannot)
}

process ComputedPaintor{
   memory params.fm_mem_req
   input :
    file(filez) from paintor_gwas
    file(ld) from ld_paintor
    file(fileannot) from paintor_fileannot
  each ncausal from 1..params.n_causal_snp
  publishDir "${params.output_dir}/paintor_$ncausal", overwrite:true, mode:'copy'
  output :
      set val(ncausal),file("${output}.results"), file("$BayesFactor") into res_paintor
  script :
    output=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_paintor_$ncausal" 
    DirPaintor=output
    annot=(params.paintor_fileannot=="" || params.paintor_annot=="") ? "" : " -annotations ${fileannot}"
    BayesFactor=output+".BayesFactor"
    """
    #mkdir -p $DirPaintor
    echo $output > input.files
    cp $filez $output
    cp $ld $output".ld"
    paint_annotation.py $fileannot $output 
    cp $fileannot $output".annotations"
    ${params.paintor_bin} -input input.files -in ./ -out ./ -Zhead Z -LDname ld -enumerate $ncausal $annot -num_samples  ${params.n_pop} -Lname $BayesFactor
    """
    
}


process ComputedCojo{
   memory params.gcta_mem_req
   cpus params.gcta_cpus_req
   input :
     set  file(bed),file(bim),file(fam) from subplink_gcta
     file(filez) from gcta_gwas
   publishDir "${params.output_dir}/cojo_gcta", overwrite:true, mode:'copy'
   output :
     file("${output}.jma.cojo")  into paintor_fileannot
     set file("${output}.cma.cojo"), file("${output}.ldr.cojo"), file("${output}.log")
   script :
    output=params.chro+"_"+params.begin_seq+"_"+params.end_seq+"_cojo"
    plk=bed.baseName
    """ 
    ${params.gcta_bin} --bfile $plk  --cojo-slct --cojo-file $filez --out $output  --cojo-p ${params.threshold_p} --thread-num ${params.gcta_cpus_req} 
    """

}
//process ComputedCojo{
//}
//-annotations "$ListHeadAnnot
//echo "$Paintor -input $DirPaintor/input.files -in $DirPaintor -out $DirPaintorOut/ -Zhead Z -LDname ld -enumerate $maxcausalpaint -annotations "$ListHeadAnnot
//cp $PWD/$rs".paintor" $DirPaintor/$rs
//cp $PWD/$rs".annot" $DirPaintor/$rs".annotations"
//cp $PWD/$rs".ld" $DirPaintor/


