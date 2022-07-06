#!/usr/bin/env nextflow
import java.nio.file.Paths

/*
 * Authors       :
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description : mtag analysis
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths
nextflow.enable.dsl = 1


// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}






def helps = [ 'help' : 'help' ]

/*params*/
allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin", "bin_mtag", "opt_mtag"]

allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2", "file_gwas"]
allowed_params+=allowed_params_head
params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}


/**/
params.boltlmm=0
/**/
params.head_pval = "P_BOLT_LMM"
params.head_freq = "A1FREQ"
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.head_beta="BETA"
params.head_se="SE"
params.head_N=""
params.head_A1="ALLELE1"
params.head_A2="ALLELE0"
params.input_dir=""
params.input_pat=""
params.list_N=""
params.mtag_mem_req="15G"
params.bin_mtag="mtag.py"
params.cut_maf = 0.001
params.opt_mtag =""
params.pheno="mtag"
params.covariates=""
params.output="output"


params.help = false

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
head_N=""

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
head_N=params.head_N
}
  
/****first step need to be format file****/
listfilegwas=params.file_gwas.split(",")
if(listfilegwas.size()<2){
       error("\n\n------\nError in your config\nFile $params.file_gwas must have at least two files \n\n---\n")
}
if(params.list_N=="" && params.input_dir=="" && params.input_pat=="" && params.head_N==""){
 error("\n\n------\nError in your config\nparams.head_N and params.list_N and params.input_dir and params.input_pat is not initialise\nnot N should be intialise\n\n---\n")


}



if(params.list_N!="" || head_N!=""){
/* Channel for each file*/
list_gwas_multi_for1=Channel.from(params.file_gwas.split(",")).flatMap{it->file(it)}
process doFormatFile{
 memory params.mtag_mem_req
 input :
   file filegwas from list_gwas_multi_for1
 output :
   file out into gwasformat
 script :
     out=filegwas.toString().replace('-','_')+".format"
     headfreq=head_freq!="" ? " --freq_header  ${head_freq} " : ""
     headN=head_freq!="" ? " --freq_N  ${head_N} " : ""
     """
     gcta_format.py --inp_asso $filegwas  --rs_header ${head_rs} --pval_header ${head_pval} $headfreq --a1_header ${head_A1} --a2_header ${head_A2} --se_header ${head_se} --beta_header ${head_beta} --chro_header ${head_chr}  --out $out --threads ${params.max_plink_cores} --bp_header ${head_bp}
     """ 
}
Head_N="N"
if( params.head_N=="")Head_N=""


}else{
list_gwas_multi_for2=Channel.from(params.file_gwas.split(",")).flatMap{it->file(it)}
bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
raw_src_ch=Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }

list_gwas_multi_for3=list_gwas_multi_for2.combine(raw_src_ch)

process doFormatFilePlk{
 memory params.mtag_mem_req
 cpus params.max_plink_cores
 input :
   set file(filegwas), file(bed), file(bim), file(fam) from list_gwas_multi_for3
 output :
   file out into gwasformat
 script :
     out=filegwas.toString().replace('-','_')+".format"
     headfreq=head_freq!="" ? " --freq_header  ${head_freq} " : ""
     headN=head_N!="" ? " --n_freq  ${head_N} " : ""
     plk=bed.baseName
     """
     gcta_format.py --inp_asso $filegwas  --rs_header ${head_rs} --pval_header ${head_pval} $headfreq --a1_header ${head_A1} --a2_header ${head_A2} --se_header ${head_se} --beta_header ${head_beta} --chro_header ${head_chr}  --out $out --threads ${params.max_plink_cores} --bp_header ${head_bp} $headN --bfile $plk  --print_pos True
     """
}
Head_N="N"
Head_freq="freq"
}
/*combine file*/
gwasformat.collect().into { listgwasform_col1; listgwasform_col2}

/*SNP bp A1 A2 freq b se N
*/

/*Mtag*/
process doMTAG{
   label 'mtag'
   memory params.mtag_mem_req
   time params.big_time
  input :
   file(listfile) from listgwasform_col1
 output :
   set val(fnames),file("${out}_trait*") into mtag_ch
 script :
   fnames = listfile.join(",")
   out='res_mtag'
   Ninfo=head_N!="" ? " --n_name ${head_N} " : " --n_value ${params.list_N}" 
   ////--info_min $info_min
   """
   ${params.bin_mtag} --sumstats $fnames --out ./$out --snp_name SNP --beta_name b --se_name se --eaf_name freq --maf_min ${params.cut_maf} --a1_name A1 --a2_name A2 --chr_name chro --p_name p  --bpos_name bp --incld_ambig_snps   ${params.opt_mtag} --force
   """
}

process RenameMtag{
   input :
     set val(listi), file(listf) from mtag_ch
   publishDir "${params.output_dir}/mtag", overwrite:true, mode:'copy'
   output :
      file("mtag_res/*") into mtag_format_ch
   script :
      listf2=listf.join(',')
      """
      mkdir -p mtag_res/
      rename_mtag.py --listi $listi --listf $listf2 --dir_out mtag_res
      """
}

mtag_format_ch_fm=mtag_format_ch.flatMap()
//SNP     CHR     BP      A1      A2      Z       N       FRQ     mtag_beta       mtag_se mtag_z  mtag_pval

process showMtag{
    memory params.other_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      file(filegwas)  from mtag_format_ch_fm
    output:
      file("${out}*")  into report_mtag
    script:
      out = filegwas.baseName.replaceAll(/_/,'-').replaceAll(/\./,'-')
      """
      general_man.py  --inp $filegwas --phenoname $out --out ${out} --chro_header CHR --pos_header BP --rs_header SNP --pval_header mtag_pval --beta_header mtag_beta --info_prog "Mtag : all file"
      """
  }

if(listfilegwas.size()>2){
list_file2by2=[]
nbfile=listfilegwas.size()
for (i = 0; i <(nbfile-1); i++){
   for( j = i+1; j<nbfile; j++){
      list_file2by2.add([i, j])
}
}

process doMTAG2by2{
   label 'mtag'
   memory params.mtag_mem_req
   time params.big_time
   input :
     file(listfile) from listgwasform_col2
     each poss from list_file2by2
     publishDir "${params.output_dir}/mtag_2by_2", overwrite:true, mode:'copy'
     output:
       file("$output"+"*")
     script :
        file1=listfile[poss[0]]
        file2=listfile[poss[1]]
        output=""+file1+"_"+file2
        """
        ${params.bin_mtag} --sumstats $file1,$file2 --out ./$output --snp_name SNP --beta_name b --se_name se --eaf_name freq --maf_min ${params.cut_maf} --a1_name A1 --a2_name A2 --chr_name chro --p_name p --bpos_name bp --incld_ambig_snps ${params.opt_mtag} --force --z_name z
        """
}
}

def getres(x) {
  def  command1 = "$x"
  def  command2 = "head -n 1"
  def proc1 = command1.execute()
  def proc2 = command2.execute()
  def proc = proc1 | proc2
  proc.waitFor()
  res ="${proc.in.text}"
  return res.trim()
}

nextflowversion =getres("nextflow -v")
if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"

report_ch = report_mtag.flatten().toList()

process doReport {
  label 'latex'
  input:
    file(reports) from report_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file("${out}.pdf")
  script:
    out = params.output+"-report"
    these_phenos     = "mtag"
    these_covariates = ""
    config = getConfig()
    images = workflow.container
    texf   = "${out}.tex"
    template "make_assoc_report.py"

}

