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



def getlistchro(listchro){
 newlistchro=[]
 for(x in listchro.split(',')) {
  splx=x.split("-")
  if(splx.size()==2){
   r1=splx[0].toInteger()
   r2=splx[1].toInteger()
   for(chro in r1..r2){
    newlistchro.add(chro.toString())
   }
  }else if(splx.size()==1){
   newlistchro.add(x)
  }else{
    logger("problem with chro argument "+x+" "+listchro)
    System.exit(0)
  }
 }
 return(newlistchro)
}


def helps = [ 'help' : 'help' ]
allowed_params = ["cut_maf", "output_dir", "pb_around_rs", "mem_req", "work_dir","mem_req","big_time", "output","nb_cpu" , "input_dir","input_pat", "file_gwas", "gwas_cat", "site_wind", "min_pval_clump", "size_win_kb"]
allowed_params_blocks = ["haploblocks", "plkref_haploblocks", "plk_othopt_haploblocks"]
allowed_params_other=["max_forks", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key","region", "AMI","maxInstances","instance-type", "instanceType", "bootStorageSize", "boot-storage-size", "max-instances", "sharedStorageMount", "shared-storage-mount", "scripts"]
allowed_params_headinfo=["head_chro_gwascat", "head_bp_gwascat", "head_pval_gwascat"]
allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2"]
allowed_params+=allowed_params_head
allowed_params+=allowed_params_other
allowed_params+=allowed_params_blocks
allowed_params+=allowed_params_headinfo
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


params.mem_req="8G"
params.big_time="1000H"

params.head_pval = "P_BOLT_LMM"
params.head_freq = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.data = ""
params.head_beta=""
params.head_se=""
params.head_A1="ALLELE1"
params.head_A2="ALLELE0"
params.other_mem_req="10GB"
params.other_cpus_req=5
params.min_pval_clump =0.001
// merge windows when take p-value
params.merge_wind=1

params.clump_r2=0.1


params.size_win_kb=250
params.size_win_kb_clump=-1
params.nb_cpu = 3

params.gwas_cat=""
// haploblocks information
params.haploblocks=""
params.plkref_haploblocks=""
params.plk_othopt_haploblocks=""

//params plink
params.plink_bin='plink'
params.genes_file=""
params.max_plink_cores = 4

params.justpheno=0

params.n_repet=10000

max_plink_cores=params.max_plink_cores 
plink_mem_req = params.plink_mem_req

plink_mem_req_max=plink_mem_req.replace('GB','000').replace('KB','').replace(' ','')
other_mem_req=params.other_mem_req
other_cpus_req=params.other_cpus_req


//params.gene_file_ftp="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"

//params gwas cat 

params.head_bp_gwascat="chromEnd"
params.head_chro_gwascat="chrom"
params.head_pval_gwascat="pValue"
params.head_af_gwascat="riskAlFreq"
params.head_riskall_gwascat="riskAllele"
params.head_pheno_gwascat="trait"
params.head_ci_gwascat="ci95"
params.head_n_gwascat="initSample"
params.head_beta_gwascat="orOrBeta"
params.head_rs_gwascat="name"
params.head_info_gwascat="pubMedID;author;trait;initSample"
//csv or tab
params.typeformat_gwascat='csv'

params.gwascat_format="USCS"

params.gwas_cat_ftp="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz"

params.pheno=""
params.file_pheno=""
params.list_chro="1-22"



//params.size_win_kb=250
//params.size_win_kb_clump=-1
size_win_kb_ld=params.size_win_kb_ld
if(size_win_kb_ld<0){
size_win_kb_ld=params.size_win_kb
}


//params.info_gwascat="DISEASE.TRAIT,REPORTED.GENE.S.,MAPPED_GENE,INITIAL.SAMPLE.SIZE"
params.threshold_pval_gwascat=1


if(params.gwas_cat==""){
process dl_gwascat_hg19{
   publishDir "${params.output_dir}/gwascat/init",  overwrite:true, mode:'copy'
   output :
      file("$out") into gwascat_init_ch
   script :
     out=params.gwas_cat_ftp.split('/')[-1]
   """
   wget -c ${params.gwas_cat_ftp}
   """

}
format="USCS"
gwascathead_chr="chrom"
gwascathead_bp="chromEnd"
infogwascat="pubMedID;author;trait;initSample"
typeformat="tab"
}else{
format="Other"
gwascat_init_ch=Channel.fromPath(params.gwas_cat,checkIfExists:true)
gwascathead_chr=params.head_chro_gwascat
gwascathead_bp=params.head_bp_gwascat
inforgwascat=params.head_info_gwascat
typeformat=params.typeformat_gwascat
infogwascat=params.head_info_gwascat
}



if(params.justpheno==1){
process get_gwascat_hg19_pheno{
   label 'R'
   publishDir "${params.output_dir}/gwascat",  overwrite:true, mode:'copy'
   input :
      file(gwascat) from gwascat_init_ch
   output :
       file("${out}*")
   script :
     out=params.output
   """
   format_gwascat_pheno.r --file $gwascat --out $out --chro_head $gwascathead_chr --bp_head $gwascathead_bp --pheno_head ${params.head_pheno_gwascat} --beta_head ${params.head_beta_gwascat} --ci_head ${params.head_ci_gwascat} --p_head ${params.head_pval_gwascat}  --n_head ${params.head_n_gwascat} --freq_head ${params.head_af_gwascat} --rs_head ${params.head_rs_gwascat} --riskall_head ${params.head_riskall_gwascat} --format $format --typeformat $typeformat
   """
}


}else{

if(params.head_beta==""){
if(params.input_dir=="" || params.input_pat==""){
println "params input_dir directory of your bedfile or/and input_pat pattern of your bedfile not define"
System.exit(-1)
}
println "beta header of your gwas file must be provide : --head_beta"
System.exit(-1)
}
if(params.head_se==""){
println "se header of your gwas file must be provide : --head_se"
System.exit(-1)
}
if(params.head_rs==""){
println "rs header of your gwas file must be provide : --head_rs"
System.exit(-1)
}
if(params.head_bp==""){
println "pos header of your gwas file must be provide : --head_bp"
System.exit(-1)
}

listchro=getlistchro(params.list_chro)
if(params.file_pheno!=''){
file_pheno_ch=file(params.file_pheno, checkIfExists:true)
}else{
file_pheno_ch=file('nofile')
}
process get_gwascat_hg19{
   label 'R'
   publishDir "${params.output_dir}/gwascat",  overwrite:true, mode:'copy'
   input :
      file(gwascat) from gwascat_init_ch
      file(filepheno) from file_pheno_ch
   output :
       file("${out}_range.bed") into gwascat_rangebed
       file("${out}.pos") into (gwascat_pos_subplk, gwascat_pos_ld2, gwascat_pos, gwascat_poswindneutre)
       file("${out}_resume.csv") into gwascat_detail_statpos
       file("${out}_all.csv") into (gwascat_all_statpos,gwascat_all_statld, gwascat_all_statclump,gwascat_all_statwind, gwascat_all_statld2,gwascat_all_statwindneutre)
       file("${out}*")
   script :
     chroparam= (params.list_chro=='') ? "" : " --chro ${listchro.join(',')}"
     phenoparam= (params.pheno=='') ? "" : " --pheno \"${params.pheno}\" "
     phenoparam= (params.file_pheno=='') ? " $phenoparam " : " --file_pheno $filepheno "
     out=params.output
   """
   format_gwascat.r --file $gwascat $chroparam $phenoparam --out $out --wind ${params.size_win_kb}  --chro_head ${params.head_chro_gwascat} --bp_head ${params.head_bp_gwascat} --pheno_head ${params.head_pheno_gwascat} --beta_head ${params.head_beta_gwascat} --ci_head ${params.head_ci_gwascat} --p_head ${params.head_pval_gwascat}  --n_head ${params.head_n_gwascat} --freq_head ${params.head_af_gwascat} --rs_head ${params.head_rs_gwascat} --riskall_head ${params.head_riskall_gwascat} --format $format

   """
}


filegwas_chrextr=Channel.fromPath(params.file_gwas,checkIfExists:true)


process extractgwas_fromgwascat{
   input : 
     file(pos) from gwascat_pos
     file(gwas) from filegwas_chrextr
   publishDir "${params.output_dir}/gwas_sub",  overwrite:true, mode:'copy'
   output :
     file("${params.output}_range.assoc") into (clump_file_ch,ld_file_ch)
     file("${params.output}_range.bed") into gwas_rangebed_subplk
     file("${params.output}_pos.init") into (pos_file_ch)
    file("${params.output}_range.init") into (range_file_ch_clump,wind_file_ch, range_file_ch_ld, range_file_ch_ld2)
     file("${params.output}*")
   script :
    wind=Math.max(size_win_kb_ld, params.size_win_kb)
    """
    extract_posgwas.py --bed $pos --gwas $gwas --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --wind $wind  --pval_gwas ${params.head_pval} --rs_gwas ${params.head_rs}  --out ${params.output} 
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


process sub_plk{
    cpus max_plink_cores
    memory plink_mem_req
  input :
      file(filegwascat) from gwascat_pos_subplk
      file(filegwas) from gwas_rangebed_subplk
      set file(bed), file(bim), file(fam) from raw_src_ch
    publishDir "${params.output_dir}/sub_plk/",  overwrite:true, mode:'copy'
    output :
       set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into (plk_ch_clump, plk_ch_ld, plk_ch_clumpstat)
    script :
    out=params.output+"_subplk"
    plkf=bed.baseName
    """
    awk '{print \$1"\t"\$2"\t"\$2"\t"\$1":"\$2}' $filegwascat > range.bed
    awk '{print \$1"\t"\$2"\t"\$2"\t"\$1":"\$2}' $filegwas >> range.bed
    plink -bfile $plkf --extract range  range.bed -out $out --make-bed --keep-allele-order --threads $max_plink_cores  --memory  $plink_mem_req_max
    """

}


process clump_aroundgwascat{
    cpus max_plink_cores
    memory plink_mem_req
   input :
      file(assocclump) from clump_file_ch
      set file(bed), file(bim), file(fam) from plk_ch_clump
   publishDir "${params.output_dir}/result/clump/tmp",  overwrite:true, mode:'copy'
   output :
      file("${out}.clumped") into clump_res_ch
      file("$out*")   
   script :
      bfile=bed.baseName
      out=params.output
      """ 
      plink -bfile $bfile  --clump $assocclump -clump-p1 $params.min_pval_clump --clump-p2 1 --clump-kb ${params.size_win_kb} --clump-r2 $params.clump_r2 -out $out --threads $max_plink_cores --memory $plink_mem_req_max
      """
}


process computedstat_pos{
   memory other_mem_req
   cpus other_cpus_req
   label 'R'
   input :
        file(assocpos) from pos_file_ch
        file(gwascat)  from gwascat_all_statpos
   publishDir "${params.output_dir}/result/exact_rep",  overwrite:true, mode:'copy'
   output :
      set file("${out}.csv"), file("${out}_cmpfrequencies.svg"), file("${out}_cmpz.svg") 
      file("$out*") 
   script :
    out=params.output+"_pos"
    af= (params.head_freq=='') ? "" : " --af_gwas ${params.head_freq} "
    """
    computestat_pos.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  $af --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out $af  --a1_gwascat ${params.head_riskall_gwascat}
    """
}


filegwas_chrextrneutre=Channel.fromPath(params.file_gwas,checkIfExists:true)
process  computedstat_windneutre{
   memory other_mem_req
   cpus other_cpus_req
   label 'R'
   input :
     file(pos) from gwascat_poswindneutre
     file(gwas) from filegwas_chrextrneutre
     file(gwascat) from  gwascat_all_statwindneutre
   publishDir "${params.output_dir}/result/wind/neutral/",  overwrite:true, mode:'copy'
   output :
      file("$output")  into random_pvalwind
   script :
    outgwas='fileassoc_clean.assoc'
    output="${params.output}_${params.size_win_kb}.pval"
    """
    extract_posgwasneutre.py --bed $pos --gwas $gwas --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --wind $params.size_win_kb --pval_gwas ${params.head_pval} --out $outgwas
    computestat_windneutre.r   --gwas $outgwas --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp}  --p_gwas $params.head_pval  --out $output --wind $params.size_win_kb  --nrep ${params.n_repet} --cpus ${params.other_cpus_req}
    """ 

}
process computedstat_win{
   memory other_mem_req
   cpus other_cpus_req
   label 'R'
   input :
        file(assocpos) from wind_file_ch
        file(gwascat)  from gwascat_all_statwind
   publishDir "${params.output_dir}/result/wind",  overwrite:true, mode:'copy'
   output :
      file("${out}*")
   script :
    out=params.output+"_wind"
    af= (params.head_freq=='') ? "" : " --af_gwas ${params.head_freq} "
    """
    computestat_wind.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  $af --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --min_pval ${params.threshold_pval_gwascat} --info_gwascat  \"$infogwascat\" --wind $params.size_win_kb --a1_gwascat ${params.head_riskall_gwascat} --merge_wind ${params.merge_wind}

    """
}



process computed_ld{
    cpus max_plink_cores
    memory plink_mem_req
  input :
      set file(bed), file(bim), file(fam) from plk_ch_ld
    publishDir "${params.output_dir}/result/ld/tmp",  overwrite:true, mode:'copy'
    output :
       file("${out}.ld") into (ld_res_ch,ld2_res_ch)
       file("$out*") 
    script :
    out=params.output+"_ld"
    plkf=bed.baseName
    """
    plink -bfile $plkf --r2  --ld-window-kb $size_win_kb_ld  --ld-window-r2 $params.clump_r2 -out $out --threads $max_plink_cores  --memory  $plink_mem_req_max  --ld-window 20000
    """

}


process computed_ld_stat{
    memory other_mem_req
    cpus other_cpus_req
    label 'R'
    input :
      file(fileld) from  ld_res_ch
      file(gwascat) from gwascat_all_statld
      file(assocpos) from range_file_ch_ld
    publishDir "${params.output_dir}/result/ld/",  overwrite:true, mode:'copy'
    output :
       file("$out*")
    script :      
      out=params.output+"_ld"
      af= (params.head_freq=='') ? "" : " --af_gwas ${params.head_freq} "
      """
      computestat_ld.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  $af --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --ld_file $fileld --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\"
      """
}

process computed_clump_stat{
    memory other_mem_req
    cpus other_cpus_req
    label 'R'
    input :
      file(fileclum) from  clump_res_ch
      file(gwascat) from gwascat_all_statclump
      file(assocpos) from range_file_ch_clump
      set file(bed), file(bim), file(fam) from plk_ch_clumpstat

    publishDir "${params.output_dir}/result/clump/",  overwrite:true, mode:'copy'
    output :
       file("$out*")
    script :
      out=params.output+"_ld"
      af= (params.head_freq=='') ? "" : " --af_gwas ${params.head_freq} "
      """
      computestat_clump.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  $af --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --clump_file $fileclum --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\" --bim $bim
      """

}
/*process computed_ld2_stat{
    memory other_mem_req
    cpus other_cpus_req
    label 'R'
    input :
      file(fileld) from  ld2_res_ch
      file(gwascatbed) from gwascat_pos_ld2
      file(gwascat) from gwascat_all_statld2
      file(assocpos) from range_file_ch_ld2
    publishDir "${params.output_dir}/result/ld2/",  overwrite:true, mode:'copy'
    output :
       file("$out")
    script :
      out_ldblock=params.output+"_ld2.tmp_pos"
      out=params.output+"_ld2"
      """
      computestat_ldv2.py  --plink_ld $fileld --pos_cat $gwascatbed --out $out
      computestat_ldv2.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --ldblock_file $out_ldblock --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\"

      """
}*/


process build_ldwind{
    memory other_mem_req
    cpus other_cpus_req
    input :
      file(fileld) from  ld2_res_ch
      file(gwascatbed) from gwascat_pos_ld2
      file(gwascat) from gwascat_all_statld2
      file(assocpos) from range_file_ch_ld2
    publishDir "${params.output_dir}/result/ldwind/tmp",  overwrite:true, mode:'copy'
    output :
      set file(gwascat), file(assocpos),file(fileld), file("${out}_ldext.out") into ldext_ch 
      set file(gwascat), file(assocpos),file(fileld), file("${out}_ldext_wind.out") into ldext_wind_ch 
      set file(gwascat), file(assocpos),file(fileld), file("${out}_ld.out") into ld_v2_ch
      set file(gwascat), file(assocpos),file(fileld), file("${out}_ld_wind.out") into ld_wind_ch 
    script :
      out=params.output+"_ldinfo"
      """
      computestat_ldbuildblock.py  --plink_ld $fileld --pos_cat $gwascatbed --out $out
      """
}

/*
process computed_ld2_stat{
   memory other_mem_req
   cpus other_cpus_req
   input :
     set file(gwascat), file(assocpos),file(fileld), file(out_ldblock) from ld_v2_ch  
    publishDir "${params.output_dir}/result/ldwind/noext",  overwrite:true, mode:'copy'
    output :
       file("$out*")
   script :
    out=params.output+"_noext"
    """
     computestat_ldv2.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --ldblock_file $out_ldblock --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\"
    """
}

process computed_ldext_stat{
   memory other_mem_req
   cpus other_cpus_req
   input :
     set file(gwascat), file(assocpos),file(fileld), file(out_ldblock) from ldext_ch
    publishDir "${params.output_dir}/result/ldwind/ext",  overwrite:true, mode:'copy'
    output :
       file("$out*")
   script :
    out=params.output+"_ext"
    """
     computestat_ldv2.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --ldblock_file $out_ldblock --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\"
    """
}

* algoritm group :
 * computed ld between positon  using `clump_r2`, and windows size of `size_win_kb`
 * defined different groups of positions
 * keep group positons where gwas catalog positon
 * extracted min p-value by group

* algoritm group extended:
 * computed ld between positon  using `clump_r2`, and windows size of `size_win_kb`
 * defined different groups of positions
 * merged groups where two positions are in two group
 * keep group positons where gwas catalog positon
 * extracted min p-value by group


*/
process computed_ldwind_stat{
 input :
      set file(gwascat), file(assocpos),file(fileld), file(out_ldwind) from ld_wind_ch
 publishDir "${params.output_dir}/result/ldwind/wind",  overwrite:true, mode:'copy'
 output:
   file("$out*")
 script :
  out_ldblock="tmpout"
  out=params.output+"_wind"
  """
  resarch_posgwas.py --gwas $assocpos --bed ${out_ldwind} --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --out ${out_ldblock}
  computestat_ldv2.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --ldblock_file $out_ldblock --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\"
  """
}

 
process computed_ldwindext_stat{
 
 input :
      set file(gwascat), file(assocpos),file(fileld), file(out_ldwind) from ldext_wind_ch
 publishDir "${params.output_dir}/result/ldwind/windext",  overwrite:true, mode:'copy'
 output:
   file("$out*")
 script :
  out_ldblock="tmpout"
  out=params.output+"_wind"
  """
  resarch_posgwas.py --gwas $assocpos --bed ${out_ldwind} --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --out ${out_ldblock}
  computestat_ldv2.r  --gwascat $gwascat --gwas $assocpos --chr_gwas ${params.head_chr} --ps_gwas ${params.head_bp} --a1_gwas ${params.head_A1} --a2_gwas ${params.head_A2}  --beta_gwas ${params.head_beta} --se_gwas ${params.head_se}  --chr_gwascat ${gwascathead_chr} --bp_gwascat ${gwascathead_bp} --p_gwas $params.head_pval --ps_gwascat $gwascathead_bp --chr_gwascat $gwascathead_chr --out $out --ldblock_file $out_ldblock --min_pvalue ${params.min_pval_clump} --min_r2  ${params.clump_r2} --info_gwascat \"$infogwascat\"
  """
}  
 

}

