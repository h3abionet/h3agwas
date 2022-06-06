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
allowed_params_head_sumstat1 = ["file_gwas_sumstat1","head_pval_sumstat1", "head_freq_sumstat1", "head_bp_sumstat1", "head_chr_sumstat1", "head_rs_sumstat1", "head_beta_sumstat1", "head_se_sumstat1", "head_A1_sumstat1", "head_A2_sumstat1", "head_n_sumstat1", "n_count1",'head_z_sumstat1']
allowed_params_head_sumstat2 = ["file_gwas_sumstat2","head_pval_sumstat2", "head_freq_sumstat2", "head_bp_sumstat2", "head_chr_sumstat2", "head_rs_sumstat2", "head_beta_sumstat2", "head_se_sumstat2", "head_A1_sumstat2", "head_A2_sumstat2", "head_n_sumstat2", 'n_count2','head_z_sumstat2']
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
if(params.head_n_sumstat1=="" & params.n_count1==""){
error(' head_n_sumstat1 not initialise and n_count1 not intitialise')
}
if(params.head_n_sumstat2=="" & params.n_count2==""){
error(' head_n_sumstat2 not initialise and n_count2 not intitialise')
}

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
    path(gwas) from filegwas_chrsum2_epos
  output:
    path(out) into sumstat2_bed
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

bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
raw_src_ch= Channel.create()
raw_src_ch2= Channel.create()

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
   .set {raw_src_ch}

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
   .set {raw_src_ch2}

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
   .set {raw_src_ch3}
plk_ch_clump  = Channel.create()
plk_ch_ld  = Channel.create()
plk_ch_clumpstat  = Channel.create()


process extractgwas_forplink{
   input : 
     path(filegwas) from filegwas_sum1
     tuple path(bed), path(bim), path(fam) from raw_src_ch2
   publishDir "${params.output_dir}/gwas_sub",  overwrite:true, mode:'copy'
   output :
     path("${out}.plink") into (file_plk_ss1_rep, file_plk_ss1, fileformat_plk1)
   script :
    out=filegwas+".assoc"
    se=(params.head_se_sumstat1=='') ? "" : " --se_header ${params.head_se_sumstat1} "
    beta=(params.head_se_sumstat1=='') ? "" : " --beta_header ${params.head_beta_sumstat1} "
    z=(params.head_z_sumstat1=='') ? "" : " --z_header ${params.head_z_sumstat1} "
    bfile=bed.baseName
    headn=(params.head_n_sumstat1=='') ? "" :  "  --n_header ${params.head_n_sumstat1} "
    ncount=(params.n_count1=='') ? "" :  "  --n ${params.n_count1} "
    freq=(params.head_freq_sumstat1=='') ? "" :  " --freq_header  ${params.head_freq_sumstat1} "
    """
    plink_format.py --inp_asso $filegwas --chro_header ${params.head_chr_sumstat1} --bp_header ${params.head_bp_sumstat1} --a1_header ${params.head_A1_sumstat1} --a2_header ${params.head_A2_sumstat1}   --pval_header ${params.head_pval_sumstat1}  --out ${out} --rs_header ${params.head_rs_sumstat1} $se $beta $z  --bfile  $bfile $headn $ncount $freq
    """
}


filegwas_sum2=Channel.fromPath(params.file_gwas_sumstat2,checkIfExists:true)
process extractgwas2_forplink{
   input : 
     path(filegwas) from filegwas_sum2
     tuple path(bed), path(bim), path(fam) from raw_src_ch3
   publishDir "${params.output_dir}/gwas_sub",  overwrite:true, mode:'copy'
   output :
     path("${out}.plink") into fileformat_plk2
   script :
    out=filegwas+".assoc"
    se=(params.head_se_sumstat2=='') ? "" : " --se_header ${params.head_se_sumstat2} "
    beta=(params.head_se_sumstat2=='') ? "" : " --beta_header ${params.head_beta_sumstat2} "
    z=(params.head_z_sumstat2=='') ? "" : " --z_header ${params.head_z_sumstat2} "
    bfile=bed.baseName
    headn=(params.head_n_sumstat2=='') ? "" :  "  --n_header ${params.head_n_sumstat2} "
    ncount=(params.n_count2=='') ? "" :  "  --n ${params.n_count2} "
    freq=(params.head_freq_sumstat2=='') ? "" :  " --freq_header  ${params.head_freq_sumstat2} "
    """
    plink_format.py --inp_asso $filegwas --chro_header ${params.head_chr_sumstat2} --bp_header ${params.head_bp_sumstat2} --a1_header ${params.head_A1_sumstat2} --a2_header ${params.head_A2_sumstat2}   --pval_header ${params.head_pval_sumstat2}  --out ${out} --rs_header ${params.head_rs_sumstat2} $se $beta $z  --bfile  $bfile $headn $ncount $freq
    """
}





process plink_sumstat1_rep{
    cpus max_plink_cores
    memory plink_mem_req
  input :
      path(filegwas) from file_plk_ss1_rep
      path(bed_pos) from sumstat2_bed 
      tuple path(bed), path(bim), path(fam) from raw_src_ch
    publishDir "${params.output_dir}/clump/",  overwrite:true, mode:'copy'
    output :
       path("${out}.clumped") into clumped_sumstat1, clumped_sumstat1_2
       path("${out}.log")
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
    path(clump) from clumped_sumstat1 
    path(sumstat1) from fileformat_plk1 
    path(sumstat2) from fileformat_plk2
  publishDir "${params.output_dir}/tmp",  overwrite:true, mode:'copy', pattern: "${out}*"
  output :
    path("$out*")
    path("${out}_all_clump.csv") into resume_all
    tuple path('sumstat_ref.metal'), path('sumstat_cmp.metal') into metalanalyse_ch
 script :
   out=params.output+"_stat" 
   gwasref="--gwas_ref $sumstat1 --gwas_ref_chr CHR --gwas_ref_bp BP --gwas_ref_a1 A1 --gwas_ref_a2 A2 --gwas_ref_beta BETA --gwas_ref_se SE --gwas_ref_p P --gwas_ref_af AF "
   gwascmp="--gwas_cmp $sumstat2 --gwas_ref_chr CHR --gwas_ref_bp BP --gwas_ref_a1 A1 --gwas_ref_a2 A2 --gwas_ref_beta BETA --gwas_ref_se SE --gwas_ref_p P --gwas_ref_af AF "
 """ 
  extract_pvalue.r $gwasref $gwascmp --file_clump $clump --out $out $n_ref $n_cmp
  """
}


process do_metal{
  label 'metaanalyse'
  input :
    tuple path(sumstat1), path(sumstat2) from metalanalyse_ch
  publishDir "${params.output_dir}/tmp/metal",  overwrite:true, mode:'copy'
  output :
     path("${out}_1.stat") into metal_res_ch
     path("${metal_config}")
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
     path(metal) from metal_res_ch
     path(fileres) from resume_all
  publishDir "${params.output_dir}/",  overwrite:true, mode:'copy'
  output : 
     path(res)
  script :
    res=params.output+'_all.csv'
    """
    merge_metalandres.r --metal $metal --resume $fileres --out $res
    """ 
}
