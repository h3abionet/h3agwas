#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths




def helps = [ 'help' : 'help' ]
allowed_params = ["cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" , "input_dir","input_pat", "file_gwas", "gwas_cat", "site_wind", "min_pval_clump", "size_win_kb","secret-key", "region", "AMI", "instanceType","bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "plink_mem_req", "plink_mem_req"]
allowed_params_other=["max_forks", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key","region", "AMI","maxInstances", "instanceType", "bootStorageSize", "boot-storage-size", "max-instances", "sharedStorageMount", "shared-storage-mount", "scripts", "max_plink_cores"]
allowed_params_head_sumstat1 = ["file_gwas_sumstat1","head_pval_sumstat1", "head_freq_sumstat1", "head_bp_sumstat1", "head_chr_sumstat1", "head_rs_sumstat1", "head_beta_sumstat1", "head_se_sumstat1", "head_A1_sumstat1", "head_A2_sumstat1", "head_n_sumstat1", "n_count1"]
allowed_params_head_sumstat2 = ["file_gwas_sumstat2","head_pval_sumstat2", "head_freq_sumstat2", "head_bp_sumstat2", "head_chr_sumstat2", "head_rs_sumstat2", "head_beta_sumstat2", "head_se_sumstat2", "head_A1_sumstat2", "head_A2_sumstat2", "head_n_sumstat2", 'n_count2']
allowed_params_clump= []


allowed_params+=allowed_params_head_sumstat1
allowed_params+=allowed_params_head_sumstat2
allowed_params+=allowed_params_other
params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}

checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}



def params_help = new LinkedHashMap(helps)
params.queue      = 'batch'
params.work_dir   = "$PWD"
params.input_listfiles = "${params.work_dir}/list_files.input"
params.output_dir = "${params.work_dir}/output"
params.output = "replication"
params.cut_maf = 0.01

params.head_n_sumstat1=""
params.n_count1=""

params.head_n_sumstat2=""
params.n_count2=""

params.mem_req="8G"
params.big_time="1000H"

params.file_gwas_sumstat1=""
params.bin_metal="metal"
params.head_pval_sumstat1 = "P_BOLT_LMM"
params.head_freq_sumstat1 = ""
params.head_bp_sumstat1 = "BP"
params.head_chr_sumstat1 = "CHR"
params.head_rs_sumstat1 = "SNP"
params.head_beta_sumstat1=""
params.head_se_sumstat1=""
params.head_A1_sumstat1="ALLELE1"
params.head_A2_sumstat1="ALLELE0"
params.head_n_sumstat1=""
params.head_ncount_sumstat1=""


params.file_gwas_sumstat2=""
params.head_pval_sumstat2 = ""
params.head_freq_sumstat2 = ""
params.head_bp_sumstat2 = ""
params.head_chr_sumstat2 = ""
params.head_rs_sumstat2 = ""
params.head_beta_sumstat2=""
params.head_se_sumstat2=""
params.head_A1_sumstat2=""
params.head_A2_sumstat2=""
params.head_n_sumstat2=""
params.head_ncount_sumstat2=""




params.input_pat=""
params.input_dir=""

params.other_cpus_req=5
params.min_pval_clump =0.001

params.clump_r2=0.5
params.clump_kb=250
params.clump_p1=0.000005
params.clump_p2=0.05


params.nb_cpu = 3

params.gwas_cat=""

//params plink
params.plink_bin='plink'
params.max_plink_cores = 4

max_plink_cores=params.max_plink_cores 
plink_mem_req = params.plink_mem_req

plink_mem_req_max=plink_mem_req.replace('GB','000').replace('KB','').replace(' ','').replace('Mb','').replace('MB','')
other_cpus_req=params.other_cpus_req



//params gwas cat 

if(params.input_dir=="" || params.input_pat==""){
println "params input_dir directory of your bedfile or/and input_pat pattern of your bedfile not define"
System.exit(-1)
}

filegwas_chrsum2_epos=Channel.fromPath(params.file_gwas_sumstat2,checkIfExists:true)
process extract_possumstat2{
  input:
    file(gwas) from filegwas_chrsum2_epos
  output:
    file(out) into sumstat2_bed
  script :
    out=gwas+".bed"
    """
    extract_chrbp_sumstat.py --inp_asso $gwas --out $out --chro_header ${params.head_chr_sumstat2}  --bp_header ${params.head_bp_sumstat2} 
    """ 
}

if(params.file_gwas_sumstat1==""){
error("\n\n------\nError in your config\nFile file_gwas_sumstat1 not defined\n\n---\n")
}
if(params.file_gwas_sumstat2==""){
error("\n\n------\nError in your config\nFile file_gwas_sumstat2 not defined\n\n---\n")
}
filegwas_sum1=Channel.fromPath(params.file_gwas_sumstat1,checkIfExists:true)

process extractgwas_forplink{
   input : 
     file(filegwas) from filegwas_sum1
   publishDir "${params.output_dir}/gwas_sub",  overwrite:true, mode:'copy'
   output :
     file("$out") into (file_plk_ss1_rep, file_plk_ss1)
   script :
    out=filegwas+".assoc"
    """
    plink_format.py --inp_asso $filegwas --chro_header ${params.head_chr_sumstat1} --bp_header ${params.head_bp_sumstat1} --a1_header ${params.head_A1_sumstat1} --a2_header ${params.head_A2_sumstat1}   --pval_header ${params.head_pval_sumstat1} --beta_header ${params.head_beta_sumstat1}  --out ${out} --rs_header ${params.head_rs_sumstat1} --se_header ${params.head_se_sumstat1}
    """
}

bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
raw_src_ch= Channel.create()

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set {raw_src_ch}

plk_ch_clump  = Channel.create()
plk_ch_ld  = Channel.create()
plk_ch_clumpstat  = Channel.create()

    /*JT : append boltlmm_assoc_ch and a]*/


process plink_sumstat1_rep{
    cpus max_plink_cores
    memory plink_mem_req
  input :
      file(filegwas) from file_plk_ss1_rep
      file(bed_pos) from sumstat2_bed 
      set file(bed), file(bim), file(fam) from raw_src_ch
    publishDir "${params.output_dir}/clump/",  overwrite:true, mode:'copy'
    output :
       file("${out}.clumped") into clumped_sumstat1, clumped_sumstat1_2
       file("${out}.log")
    script :
    out=params.output+"_subplk"
    plkf=bed.baseName
    """
    plink -bfile $plkf --extract range ${bed_pos}  -out $out  --keep-allele-order --threads $max_plink_cores  --memory  $plink_mem_req_max  --clump $filegwas --clump-p1 ${params.clump_p1} --clump-p2 ${params.clump_p2} --clump-r2 ${params.clump_r2} --clump-kb ${params.clump_kb}
    """
}

filegwas_sum1_stat=Channel.fromPath(params.file_gwas_sumstat1,checkIfExists:true)
filegwas_sum2_stat=Channel.fromPath(params.file_gwas_sumstat2,checkIfExists:true)


process computed_stat{
  label 'R'
  input :
    file(clump) from clumped_sumstat1 
    file(sumstat1) from filegwas_sum1_stat
    file(sumstat2) from filegwas_sum2_stat
  publishDir "${params.output_dir}/tmp",  overwrite:true, mode:'copy', pattern: "${out}*"
  output :
    file("$out*")
    file("${out}_all_clump.csv") into resume_all
    set file('sumstat_ref.metal'), file('sumstat_cmp.metal') into metalanalyse_ch
 script :
   out=params.output+"_stat" 
   n_ref=(params.head_n_sumstat1!="") ? " --gwas_ref_n ${params.head_n_sumstat1} " :  " --gwas_ref_ncount ${params.head_ncount_sumstat1}"
   n_cmp=(params.head_n_sumstat2!="") ? " --gwas_cmp_n ${params.head_n_sumstat2} " :  " --gwas_cmp_ncount ${params.head_ncount_sumstat2}"
 """ 
#
 extract_pvalue.r \
 --gwas_ref $sumstat1 --gwas_ref_chr ${params.head_chr_sumstat1} --gwas_ref_bp ${params.head_bp_sumstat1} --gwas_ref_a1 ${params.head_A1_sumstat1} --gwas_ref_a2 ${params.head_A2_sumstat1} --gwas_ref_beta ${params.head_beta_sumstat1} --gwas_ref_se ${params.head_se_sumstat1} --gwas_ref_p ${params.head_pval_sumstat1} --gwas_ref_af ${params.head_freq_sumstat1} \
 --gwas_cmp $sumstat2 --gwas_cmp_chr ${params.head_chr_sumstat2} --gwas_cmp_bp ${params.head_bp_sumstat2} --gwas_cmp_a1 ${params.head_A1_sumstat2} --gwas_cmp_a2 ${params.head_A2_sumstat2} --gwas_cmp_beta ${params.head_beta_sumstat2} --gwas_cmp_se ${params.head_se_sumstat2} --gwas_cmp_p ${params.head_pval_sumstat2} --gwas_cmp_af ${params.head_freq_sumstat2} \
 --file_clump $clump --out $out $n_ref $n_cmp
  """
}


process do_metal{
  label 'metaanalyse'
  input :
    set file(sumstat1), file(sumstat2) from metalanalyse_ch
  publishDir "${params.output_dir}/tmp/metal",  overwrite:true, mode:'copy'
  output :
     file("${out}1.stat") into metal_res_ch
     file("${metal_config}")
  script :
      metal_config="metal_config.config"
      out=params.output+'_metal'
      """
      echo $sumstat1 > fileListe
      echo $sumstat2 >> fileListe
      ma_get_configmetal.py --filelist fileListe  --output_configmetal $metal_config  --out_file_metal $out --genomic_control F
      ${params.bin_metal} $metal_config
      """
}

process merge_res{
  label 'R'
  input  :
     file(metal) from metal_res_ch
     file(fileres) from resume_all
  publishDir "${params.output_dir}/",  overwrite:true, mode:'copy'
  output : 
     file(res)
  script :
    res=params.output+'_all.csv'
    """
    merge_metalandres.r --metal $metal --resume $fileres --out $res
    """ 
}
