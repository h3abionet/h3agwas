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

nextflow.enable.dsl = 1

if (!workflow.resume) {
    def dir = new File(params.output_dir)
    if (dir.exists() && dir.directory && (!(dir.list() as List).empty)) {
       println "\n\n============================================"
       println "Unless you are doing a -resume, the output directory should be empty"
       println "We do not want to overwrite something valuable in "+params.output_dir
       println "Either clean your output directory or check if you meant to do a -resume"
       System.exit(-1)
    }
}


def strmem(val){
 return val as nextflow.util.MemoryUnit
}

/*definition*/
def errormess(message,exitn=-1){
    if(message=="")return(0)
    println(message)
    System.exit(exitn)
}
def checkparams(param, namesparam, type, min=null, max=null, possibleval=null, notpossibleval=null) {
  messageerror=""
  if(param==null){
    messageerror+="error :--"+namesparam+" is null " 
  } else {
    if(!(param.getClass() in type)){
   messageerror+="error :--"+namesparam+" must be a "+ type
     if(params.getClass()==Boolean)messageerror+=", but no parameters given"
     else messageerror+=" but type is "+param.getClass()+" value "+ param
   }else{
   if(min && param<min)messageerror+="\nerror : --"+namesparam+" < min value :"+param +" < "+min
   if(max && param>max)messageerror+="\nerror : --"+namesparam +"> maxvalue :" + param+" > "+max
   if(possibleval && !(param in possibleval))messageerror+="\nerro : --"+namesparam +" must be one the value :"+possibleval.join(',')
   }
   }
   errormess(messageerror,3)
}


def checkmultiparam(params, listparams, type, min=null, max=null, possibleval=null, notpossibleval=null){
 messageerror=""
 for(param in listparams){
   if(params.containsKey(param)){
     checkparams(params[param], param, type, min=min, max=max, possibleval=possibleval, notpossibleval=notpossibleval)
   }else{
     messageerror+="param :"+param+" not initialize\n" 
   }
 }
 errormess(messageerror, 2)
}



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
    System.exit(-1)
  }
 }
 return(newlistchro)
}
//---- General definitions --------------------------------------------------//

import java.nio.file.Paths

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}
nullfile = [false,"False","false", "FALSE",0,"","0","null",null]
def checkColumnHeader(fname, columns) {
  if (workflow.profile == "awsbatch") return;
  if (fname.toString().contains("s3://")) return;
  if (fname.toString().contains("az://")) return;
  if (nullfile.contains(fname)) return;
  new File(fname).withReader { line = it.readLine().tokenize() }
  problem = false;
  columns.each { col ->
    if (! line.contains(col) & col!='') {
      println "The file <$fname> does not contain the column <$col>";
      problem=true;
    }
    if (problem)
      System.exit(-1)
  }
}





def helps = [ 'help' : 'help' ]

allowed_params_input = ["input_dir","input_pat","output","output_dir","plink_mem_req", "work_dir", "scripts",  "accessKey", "access-key", "secretKey", "secret-key", "region",  "big_time",  "rs_list", 'list_phenogc', 'cojo_slct_other', "paintor_fileannot", "paintor_listfileannot", "caviarbf_avalue", "gwas_cat", "genes_file", "genes_file_ftp", "list_phenogc", "file_phenogc", "headgc_chr", "headgc_bp", "headgc_bp", "genes_file","genes_file_ftp", "list_chro",  'modelsearch_caviarbf_bin', "AMI", "instanceType", "instance-type", "bootStorageSize","maxInstances", "max-instances", "sharedStorageMount", "shared-storage-mount",'queue', "data", 'pheno', 'covariates']
allowed_params=allowed_params_input
allowed_params_bin=["finemap_bin", "paintor_bin","plink_bin", "caviarbf_bin", "gcta_bin", "gwas_cat_ftp"]
allowed_params+=allowed_params_bin
allowed_params_cores=["plink_cpus_req", "gcta_cpus_req", "fm_cpus_req",'max_plink_cores']
allowed_params+=allowed_params_cores
allowed_params_intother=["max_forks", "n_pop", "n_causal_snp", 'cojo_slct', "size_wind_kb"]
allowed_params+=allowed_params_intother
allowed_params_bolother=["used_pval_z"]
allowed_params+=allowed_params_bolother
allowed_params_float=["cut_maf", "threshold_p", "prob_cred_set", 'threshold_p2', 'clump_r2']
allowed_params+=allowed_params_float
allowed_params_memory=["gcta_mem_req" , "plink_mem_req", "other_mem_req", "fm_mem_req","caviar_mem_req", "paintor_mem_req", "boot-storage-size"]
allowed_params+=allowed_params_memory
params_filegwas=[ "file_gwas", "head_beta", "head_se", "head_A1", "head_A2", "head_freq", "head_chr", "head_bp", "head_rs", "head_pval", "head_n"]
allowed_params+=params_filegwas




def params_help = new LinkedHashMap(helps)


filescript=file(workflow.scriptFile)
projectdir="${filescript.getParent()}"
dummy_dir="${projectdir}/../qc/input"


params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.genes_file=""
params.genes_file_ftp="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
params.output="finemap"

params.gcta_bin="gcta64"

// paramater
params.n_pop=0

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
params.head_n=""
params.data=""
params.cut_maf=0.01
params.used_pval_z=0
params.headgc_chr=""
params.headgc_bp=""
params.gwas_cat = ""
params.file_phenogc = ""

params.file_gwas=""

params.prob_cred_set=0.95
params.covariates=""

params.other_mem_req="20GB"

// gcta parameters
params.gcta_mem_req="15GB"
params.gcta_cpus_req = 1
params.plink_cpus_req=5

params.fm_cpus_req = 5
params.cojo_slct=1
params.cojo_slct_other=""

params.caviar_mem_req="40GB"
params.paintor_mem_req="20GB"
params.fm_mem_req = "20G"
params.plink_mem_req="6GB"


params.big_time='100h'

params.threshold_p=5*10**-8
params.n_causal_snp=1
params.caviarbf_avalue="0.1,0.2,0.4"
params.paintor_fileannot=""
params.paintor_listfileannot=""
params.threshold_p2=0.5
params.clump_r2=0.5
params.size_wind_kb=100
//params.paintor_annot=""

params.gwas_cat_ftp="http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz"
params.list_chro="1-22"
params.list_phenogc=""
params.pheno=""




params.finemap_bin="finemap"
params.caviarbf_bin="caviarbf"
params.modelsearch_caviarbf_bin="model_search"
params.paintor_bin="PAINTOR"
params.plink_bin="plink"

listchro=getlistchro(params.list_chro)
if(params.gwas_cat==""){
println('gwas_cat : gwas catalog option not initialise, will be downloaded')
if(params.file_phenogc=="")phenogc_ch=channel.fromPath("${dummy_dir}/01")
else phenogc_ch=channel.fromPath(params.file_phenogc)

process GwasCatDl{
    label 'R'
    publishDir "${params.output_dir}/gwascat",  mode:'copy'
    input :
      file(phenogc) from  phenogc_ch
    output :
       file("${out}_all.csv") into gwascat_ch
       file("${out}*")
    script :
      phenol= (params.list_phenogc=="") ? "" : "  --pheno '${params.list_phenogc}' "
      phenofile= (params.file_phenogc=="") ? "" : "  --file_pheno  $phenogc "
      out="gwascat_format"
      """
      wget -c ${params.gwas_cat_ftp} --no-check-certificate
      format_gwascat.r --file `basename ${params.gwas_cat_ftp}` $phenol --out $out  --chro ${listchro.join(',')} $phenofile
      """ 
}
headgc_chr="chrom"
headgc_bp="chromEnd"
}else{
gwascat_ch=Channel.fromPath(params.gwas_cat, checkIfExists:true)
headgc_chr=params.headgc_chr
headgc_bp=params.headgc_bp

}

if(params.file_gwas==""){
    println "\nfile_gwas not initialised\n";
    System.exit(-1)
}


params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}

checkmultiparam(params,allowed_params_input, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_memory, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
//checkmultiparam(params,allowed_params_bin, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_cores, java.lang.Integer, min=1, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_intother, java.lang.Integer, min=0, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_bolother, java.lang.Integer, min=0, max=null, possibleval=[0,1], notpossibleval=null)
checkmultiparam(params,allowed_params_float, [java.lang.Double,java.lang.Float, java.lang.Integer, java.math.BigDecimal], min=0, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,params_filegwas, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)


bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")

raw_src_ch= Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


//gwas_extract_plk=Channel.create()
//plink_subplk=Channel.create()
//gwas_plk_clump=Channel.create()
if(params.data){
data_ch=channel.fromPath(params.data, checkIfExists:true)
process extract_inddata{
  label 'R'
  input :
    path(data) from data_ch
  output :
     path(newkeep) into file_keep
  script :
       newkeep=data+"_keep"
       pheno=(params.pheno=="") ? "" : " --pheno ${params.pheno}"
       cov=(params.covariates=="") ? "" : " --cov ${params.covariates} "
       """
       extract_indplink.r --data $data $pheno --out $newkeep $cov
       """
}

}else {
file_keep=channel.fromPath("$dummy_dir/00")
}
process clean_plink{
   memory params.plink_mem_req
   cpus params.max_plink_cores
   input :
     tuple path(bed), path(bim), path(fam)  from raw_src_ch
     path(filekeep) from file_keep
   output :
     tuple path("${newbed}.bed"),path("${newbed}.bim"),path("${newbed}.fam") into gwas_extract_plk, plink_subplk,gwas_plk_clump
   script :
     bfile=bed.baseName
     newbed=bfile+"_clean"
     maxmaf=1 - params.cut_maf
     keep=(params.data=="") ? "" : " --keep $filekeep "
     """ 
     plink -bfile $bfile --maf ${params.cut_maf}  -out  $newbed -make-bed --threads ${params.max_plink_cores} $keep
     """
}
//raw_src_ch.separate( gwas_extract_plk, plink_subplk,gwas_plk_clump) { a -> [ a, a,a] }


gwas_file=Channel.fromPath(params.file_gwas,checkIfExists:true)
gwas_file_clump=Channel.fromPath(params.file_gwas,checkIfExists:true)
// plink 
checkColumnHeader(params.file_gwas, [params.head_beta, params.head_se, params.head_A1,params.head_A2, params.head_freq, params.head_chr, params.head_bp, params.head_rs, params.head_pval, params.head_n])
process clump_data{
 memory params.plink_mem_req
 cpus params.max_plink_cores
 input :
     set file(bed),file(bim),file(fam) from gwas_plk_clump
     file(gwasfile) from gwas_file_clump
 publishDir "${params.output_dir}/clump/",  mode:'copy'
 output :
    file("${output}.clumped") into (file_clump,file_clump2)
 script :
   output="clump_output"
   plkfile=bed.baseName
   plink_mem_req_max=params.plink_mem_req.replace('GB','000').replace('KB','').replace(' ','').replace('MB','').replace('Mb','')
   afhead=(params.head_freq=='') ? "" : " --freq_header ${params.head_freq} "  
   """
   formatsumstat_inplink.py --inp_asso $gwasfile --chro_header  ${params.head_chr} --bp_header ${params.head_bp} --a1_header ${params.head_A1} --a2_header ${params.head_A2}  --pval_header ${params.head_pval} --beta_header ${params.head_beta}  --out $output --rs_header ${params.head_rs} --se_header  ${params.head_se} --bim $bim --maf ${params.cut_maf} $afhead
 ${params.plink_bin} -bfile $plkfile  -out $output  --threads ${params.max_plink_cores}   --clump $output --clump-p1 ${params.threshold_p} --clump-p2 ${params.threshold_p2} --clump-r2 ${params.clump_r2} --clump-kb ${params.size_wind_kb} --maf ${params.cut_maf} --memory $plink_mem_req_max
  """
}

process extract_sigpos_gwas{
  memory params.plink_mem_req
  memory params.max_plink_cores
  input :
     path(filegwas) from gwas_file
     tuple path(bed),path(bim),path(fam) from gwas_extract_plk
     path(clump) from file_clump2
  output :
    file("*.gcta") into gcta_gwas_i
    file("*_finemap.z") into  (finemap_gwas_cond_i, finemap_gwas_sss_i)
    file("*_caviar.z") into caviarbf_gwas_i
    file("*.paintor") into paintor_gwas_i
    file("*.range") into range_plink_i
    file("*.all") into data_i_i
    file("*.pos") into paintor_gwas_annot_i
    env n into nval_ch, nval_ch_fmss, nval_ch_fm, nval_ch_caviar, nval_ch_paintor
  script :
    freq= (params.head_freq=="") ? "":" --freq_header ${params.head_freq} "
    nheader= (params.head_n=="") ? "":" --n_header ${params.head_n}"
    nvalue= (params.n_pop=="") ? "":" --n ${params.n_pop}"
    bfile=bed.baseName
    around=params.size_wind_kb*1000
    out= params.output
    """
    fine_extract_sig_mw.py --inp_resgwas $filegwas --chro_header ${params.head_chr} --pos_header ${params.head_bp} --beta_header ${params.head_beta} --se_header ${params.head_se} --a1_header ${params.head_A1} --a2_header ${params.head_A2} $freq  --bfile $bfile --rs_header ${params.head_rs} --out_head $out --p_header ${params.head_pval}  $nvalue --min_pval ${params.threshold_p} $nheader --z_pval ${params.used_pval_z} --clump $clump --around $around  --maf ${params.cut_maf} --threads ${params.max_plink_cores}
    n=`cat n.out`
    """
}

///*
process extract_sigpos{
  input :
   file(clump) from file_clump
  output :
     stdout into  postonalyse2
  script:
  """
   sed '1d' $clump| awk '{if(\$1!="")print \$1"_"\$4}' 
  """
}

range_plink=range_plink_i.flatMap{it}
range_plink_ch=range_plink.combine(plink_subplk)

process SubPlink{
  memory params.plink_mem_req
  cpus params.max_plink_cores
  input :
     set file(range), file(bed),file(bim),file(fam) from range_plink_ch
  output :
     set val(pos),file("${out}.bed"),file("${out}.bim"),file("${out}.fam") into (subplink_ld, subplink_gcta)
  script : 
     pos=range.baseName
     plink_mem_req_max=params.plink_mem_req.replace('GB','000').replace('KB','').replace(' ','').replace('MB','').replace('Mb','')
     plk=bed.baseName
     out=pos+'_sub'
     """
     ${params.plink_bin} -bfile $plk  --keep-allele-order --extract  range  $range --make-bed -out  $out  # --threads ${params.max_plink_cores}  --memory  $plink_mem_req_max
     """
}

process ComputedLd{
   memory params.plink_mem_req
   cpus params.max_plink_cores
   input : 
      set val(pos),file(bed),file(bim),file(fam) from subplink_ld
  output :
       set val(pos),file("$outld") into (ld_fmcond, ld_fmsss,ld_caviarbf, ld_paintor)
   script :
  plink_mem_req_max=params.plink_mem_req.replace('GB','000').replace('KB','').replace(' ','').replace('MB','').replace('Mb','')
    outld=pos+".ld"
    plk=bed.baseName
    """
    ${params.plink_bin} --r2 square0 yes-really -bfile $plk -out "tmp" --threads ${params.max_plink_cores}  --memory  $plink_mem_req_max
    sed 's/\\t/ /g' tmp.ld | sed 's/nan/0/g' > $outld
    """
}
finemap_gwas_cond_2=finemap_gwas_cond_i.flatMap{it}
process format_ldfmcond{
   input :
      file(finem) from finemap_gwas_cond_2 
   output :
      tuple val(pos), path(finem) into finemap_gwas_cond
   script :
       spl=finem.baseName.split('_')
       pos=spl[0]+'_'+spl[1]
       """
       echo $pos
       """
}

ld_fmcond_group=ld_fmcond.join(finemap_gwas_cond).combine(nval_ch_fm)

process ComputedFineMapCond{
  label 'finemapping'
  memory params.fm_mem_req
  input :
    set val(pos),file(ld),file(filez), val(n) from ld_fmcond_group
  publishDir "${params.output_dir}/$pos/fm_cond",  mode:'copy'
  output :
    set val(pos), file("${out}.snp"),file("${out}.cred") into res_fmcond
   set file("${out}.config"), file("${out}.cred"), file("${out}.log_cond")
  script:
  fileconfig="config"
  out=pos+"_cond" 
  """ 
  echo "z;ld;snp;config;cred;log;n_samples" > $fileconfig
  echo "$filez;$ld;${out}.snp;${out}.config;${out}.cred;${out}.log;${n}" >> $fileconfig
  ${params.finemap_bin} --cond --in-files $fileconfig   --log --cond-pvalue ${params.threshold_p}  --n-causal-snps ${params.n_causal_snp}  --prob-cred-set ${params.prob_cred_set} 
  """
}

finemap_gwas_sss_2=finemap_gwas_sss_i.flatMap{it}
process format_ldfmsss{
   input :
      file(finem) from finemap_gwas_sss_2
   output :
      tuple val(pos), path(finem) into finemap_gwas_sss
   script :
       spl=finem.baseName.split('_')
       pos=spl[0]+'_'+spl[1]
       """
       echo $pos
       """
}

ld_fmss_group=ld_fmsss.join(finemap_gwas_sss).combine(nval_ch_fmss)


process ComputedFineMapSSS{
  label 'finemapping'
  memory params.fm_mem_req
  cpus params.fm_cpus_req
  input :
     set val(pos),file(ld),file(filez),val(n) from ld_fmss_group
  publishDir "${params.output_dir}/$pos/fm_sss",  mode:'copy'
  output :
    set val(pos),file("${out}.snp"), file("${out}.cred${params.n_causal_snp}") into res_fmsss
    set file("${out}.config"), file("${out}.cred${params.n_causal_snp}"), file("${out}.log_sss")
  script:
  fileconfig="config"
  out=pos+"_sss"
  """
  echo "z;ld;snp;config;cred;log;n_samples" > $fileconfig
  echo "$filez;$ld;${out}.snp;${out}.config;${out}.cred;${out}.log;$n" >> $fileconfig
  ${params.finemap_bin} --sss --in-files $fileconfig  --n-threads ${params.fm_cpus_req}  --log --n-causal-snps ${params.n_causal_snp} --prob-cred-set ${params.prob_cred_set}
  """
}

caviarbf_gwas_2=caviarbf_gwas_i.flatMap{it}
process format_caviarbf{
   input :
      file(finem) from caviarbf_gwas_2
   output :
      tuple val(pos), path(finem) into caviarbf_gwas
   script :
       spl=finem.baseName.split('_')
       pos=spl[0]+'_'+spl[1]
       """
       echo $pos
       """
}


ld_caviarbf_group=ld_caviarbf.join(caviarbf_gwas).combine(nval_ch_caviar)
process ComputedCaviarBF{
  memory params.caviar_mem_req
  label 'finemapping'
  input :
     set val(pos),file(ld),file(filez), val(n) from ld_caviarbf_group
  publishDir "${params.output_dir}/$pos/caviarbf",  mode:'copy'
  output :
   set val(pos), file("${output}.marginal") into res_caviarbf
   set file("$output"), file("${output}.statistics")
  script :
   output=pos+"_caviarbf"
   """
   ${params.caviarbf_bin} -z ${filez} -r $ld  -t 0 -a ${params.caviarbf_avalue} -c ${params.n_causal_snp} -o ${output} -n $n
   nb=`cat ${filez}|wc -l `
   ${params.modelsearch_caviarbf_bin} -i $output -p 0 -o $output -m \$nb 2> ${output}_modelsearch.log
   """
}

paintor_gwas_annot_2=paintor_gwas_annot_i.flatMap{it}
process paintor_format{
   input :
      file(finem) from paintor_gwas_annot_2
   output :
      tuple val(pos), path(finem) into paintor_gwas_annot
      val(pos) into (list_pos_paint_ch1, list_pos_paint_ch2)
   script :
       spl=finem.baseName.split('_')
       pos=spl[0]+'_'+spl[1]
       """
       echo $pos
       """
}


NCausalSnp=Channel.from(1..params.n_causal_snp)
baliseannotpaint=0
if(params.paintor_bin!="0" & params.paintor_bin!=0 & params.paintor_bin!=""){
 if(params.paintor_fileannot!=""){
  paintor_fileannot=pos_tonalyse_ch_1.combine(Channel.fromPath("${dummy_dir}/01", checkIfExists:true))
  paintor_fileannotplot=pos_tonalyse_ch_2.combine(Channel.fromPath("${dummy_dir}/02"))
  baliseannotpaint=1
 }else{
  if(params.paintor_listfileannot!=""){
   baliseannotpaint=1
   paintor_gwas_annot2=paintor_gwas_annot.combine(Channel.fromPath(params.paintor_listfileannot))
   process paintor_selectannot{
    input :
     set val(pos),file(list_loc), file(listinfo) from paintor_gwas_annot2
    publishDir "${params.output_dir}/$pos/paintor/annot", mode:'copy'
    output :
     set val(pos),file(out) into (paintor_fileannot, paintor_fileannotplot, paintor_fileannot2)
    script :
     outtmp="tmp.res"
     out="annotationinfo"
     """
     head -1 $list_loc > $outtmp
     sed '1d' $list_loc |awk '{print "chr"\$0}' >> $outtmp
    annotate_locus_paint.py --input $listinfo  --locus $outtmp --out $out --chr chromosome --pos position
    """
  }
  process paintor_extractannotname{
    input :
       set val(pos),file(fileannot) from paintor_fileannot2
    output :
       set val(pos),stdout into annotname
    """
    head -1 $fileannot | sed 's/ /,/g' 
    """ 
  }
 } else{
  println 'no file annot for paintor'
  postonalyse2.into{postonanalyse_tmp1; postonanalyse_tmp2;postonanalyse_tmp3}
  pos_tonalyse_ch_1=postonanalyse_tmp1.flatMap { list_str -> list_str.split() }
  paintor_fileannot=pos_tonalyse_ch_1.combine(Channel.fromPath("${dummy_dir}/00", checkIfExists:true))
  pos_tonalyse_ch_2=postonanalyse_tmp2.flatMap { list_str -> list_str.split() }
  paintor_fileannotplot=pos_tonalyse_ch_2.combine(Channel.fromPath("${dummy_dir}/01"))
  pos_tonalyse_ch_3=postonanalyse_tmp3.flatMap { list_str -> list_str.split() }
  annotname=pos_tonalyse_ch_3.combine(Channel.value("N"))
  }
 }
 
paintor_gwas_2=paintor_gwas_i.flatMap{it}
 process paintor_format2{
   input :
      file(finem) from paintor_gwas_2
   output :
      tuple val(pos), path(finem) into paintor_gwas
   script :
       spl=finem.baseName.split('_')
       pos=spl[0]+'_'+spl[1]
       """
       echo $pos
       """
 }
 //paintor_gwas=ld_paintor.join(paintor_gwas)

 ld_paintor_group=ld_paintor.join(paintor_gwas).join(paintor_fileannot).join(annotname).combine(nval_ch_paintor)// .combine(annotname).combine(paintor_fileannot).combine(nval_ch_paintor)
 process ComputedPaintor{
   label 'finemapping'
   memory params.paintor_mem_req
   input :
    set val(pos),file(ld),file(filez), file(fileannot), val(annot_name),val(n) from ld_paintor_group
  publishDir "${params.output_dir}/$pos/paintor/",  mode:'copy'
  output :
      set val(pos),file("${output}.results") into res_paintor_ch
      set val(pos),file(FileInfo) into infores_paintor_ch
      file("${output}*")
      set val(pos), file(BayesFactor) into res_paintor_ch_bf
  script :
    ncausal=params.n_causal_snp
    output=pos+"_paintor_$ncausal" 
    DirPaintor=output
    annot=(baliseannotpaint==0) ? "" : " -Gname ${output}_an  -annotations ${annot_name}"
    BayesFactor=output+".BayesFactor"
    FileInfo=output+".info"
    Info="$ncausal;${output}.results;$BayesFactor"
    """
    echo "$Info" > $FileInfo
    echo $output > input.files
    cp $filez $output
    cp $ld $output".ld"
    if [ "$fileannot" == "00" ]
    then
    paint_annotation.py $fileannot $output  $output".annotations"
    else 
    cp $fileannot $output".annotations"
    fi
    ${params.paintor_bin} -input input.files -in ./ -out ./ -Zhead Z -LDname ld -enumerate $ncausal -num_samples  $n -Lname $BayesFactor $annot
    """
 }
}else{
 ld_paintor_group=ld_paintor.combine(channel.from("00"))
 process ComputedPaintor_null{
  input :
   set val(pos), file(ld),file(filez) from ld_paintor_group
  output :
     set val(pos),file("${output}.results") into res_paintor_ch
     set val(pos),file(FileInfo) into infores_paintor_ch
     set val(pos), file(BayesFactor) into res_paintor_ch_bf
  script :
    output=pos+"_paintor_null"
    FileInfo=output+".info"
    BayesFactor=output+".BayesFactor"
     """
     touch $output".results"
     touch $BayesFactor
     touch $FileInfo
     """
 }
paintor_fileannot=list_pos_paint_ch1.combine(channel.fromPath("${dummy_dir}/02"))
paintor_fileannotplot=list_pos_paint_ch2.combine(channel.fromPath("${dummy_dir}/03"))
annotname=Channel.from("N")

}


 gcta_gwas_2=gcta_gwas_i.flatMap{it}
 process gcta_format{
   input :
      file(finem) from gcta_gwas_2
   output :
      tuple val(pos), path(finem) into gcta_gwas
   script :
       spl=finem.baseName.split('_')
       pos=spl[0]+'_'+spl[1]
       """
       echo $pos
       """
 }
 //

gcta_gwas_join=gcta_gwas.join(subplink_gcta)
process ComputedCojo{
   label 'gcta'
   memory params.gcta_mem_req
   cpus params.gcta_cpus_req
   input :
     set val(pos),file(filez), file(bed),file(bim),file(fam) from gcta_gwas_join
   publishDir "${params.output_dir}/$pos/cojo_gcta",  mode:'copy'
   output :
     set val(pos), file("${output}.jma.cojo")  into res_cojo
     set file("${output}.cma.cojo"), file("${output}.ldr.cojo"), file("${output}.log")
   script :
    output=pos+"_cojo"
    plk=bed.baseName
    """ 
    ${params.gcta_bin} --bfile $plk  --cojo-slct --cojo-file $filez --out $output  --cojo-p ${params.threshold_p} --thread-num ${params.gcta_cpus_req}  --diff-freq 0.49
    if [ ! -f $output".jma.cojo" ]
    then 
    echo -e  "Chr	SNP	bp	refA	freq	b	se	p	n	freq_geno	bJ	bJ_se	pJ	LD_r"	> $output".jma.cojo"
    touch $output".cma.cojo" $output".ldr.cojo"
    fi
    """
}
if(params.genes_file==""){
process GetGenesInfo{
   memory { strmem(params.other_mem_req) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 10

   output :
      file(out) into genes_file_ch
   publishDir "${params.output_dir}/data/", mode:'copy'
   script :
     out="gencode.v19.genes"
     """
     wget -c ${params.genes_file_ftp} --no-check-certificate
     zcat `basename ${params.genes_file_ftp}` > file_genes
     change_genes_gencode.py file_genes
     """
}
}else{
genes_file_ch=Channel.fromPath(params.genes_file)
}

 data_i2=data_i_i.flatMap{it}
 process datai_format{
   input :
      file(finem) from data_i2
   output : 
      tuple val(pos), path(finem) into data_i
   script :
       spl=finem.baseName.split('_')
       pos=spl[0]+'_'+spl[1]
       """
       echo $pos
       """
 }

mergeall=res_paintor_ch.join(infores_paintor_ch).join(res_paintor_ch_bf).join(paintor_fileannotplot).join(res_cojo).join(res_caviarbf).join(res_fmsss).join(res_fmcond).join(data_i).combine(genes_file_ch).combine(gwascat_ch)
process MergeResult{
    label 'R'
    memory params.other_mem_req
    input :
      set val(pos), file(paintori),file(infopaintor), file(paintor_bf), file(pfileannoti), file(cojo), file(caviarbf), file(fmsss), file(fmssscred), file(fmcond), file(fmcondcred), file(datai), file(genes), file(gwascat) from mergeall
   publishDir "${params.output_dir}/$pos/",   pattern:"${out}*", mode:'copy'
   publishDir "${params.output_dir}/$pos/plot",    mode:'copy'
    output :
       set file("${out}.pdf"), file("${out}.all.out"), file("${out}.all.out")
       set val(pos), file(paintori),file(infopaintor), file(paintor_bf), file(pfileannoti), file(cojo), file(caviarbf), file(fmsss),file(fmssscred), file(fmcond),file(fmcondcred), file(datai), file(genes), file(gwascat), file("run.bash"), file("infopaintor"), file("infofinemap"), file("infofinemapcred")
    script :
      out=params.output+'_'+pos
      infopaint=infopaintor.join(" ")
      pfileannot= (baliseannotpaint=="0" | params.paintor_bin==0 | params.paintor_bin=="" ) ? "":" --paintor_fileannot $pfileannoti "
      paintor = (params.paintor_bin==0 | params.paintor_bin=="") ? "" : "--listpaintor  infopaintor"
      
      """
       cat $infopaintor > infopaintor
       echo "sss $fmsss" > infofinemap 
       echo "sss $fmssscred" > infofinemapcred
       echo "cond $fmcond" >> infofinemap 
       echo "cond $fmcondcred" >> infofinemapcred
       merge_finemapping_v2.r --out $out $paintor  --cojo  $cojo --datai  $datai --caviarbf $caviarbf --list_genes $genes  --gwascat $gwascat --headbp_gc ${headgc_bp} --headchr_gc ${headgc_chr}  --listfinemap infofinemap  $pfileannot --listfinemap_cred infofinemapcred
       echo "merge_finemapping_v2.r --out $out --listpaintor  infopaintor  --cojo  $cojo --datai  $datai --caviarbf $caviarbf --list_genes $genes  --gwascat $gwascat --headbp_gc ${headgc_bp} --headchr_gc ${headgc_chr}  --listfinemap infofinemap  $pfileannot --listfinemap_cred infofinemapcred" > run.bash
      """
}

