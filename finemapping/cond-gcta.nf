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
 * Description : pipeline to do a Conditional and joint multiple-SNP analysis of GWAS
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths
/*definition*/
nextflow.enable.dsl = 1
def errormess(message,exitn=0){
    if(message=="")return(0)
    println(message)
    System.exit(exitn)
}
def strmem(val){
 return val as nextflow.util.MemoryUnit
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
    errormess(messageerror,-1)
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
 errormess(messageerror, -1)
}

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







def helps = [ 'help' : 'help' ]


allowed_params_input = ["input_dir","input_pat","output","output_dir","plink_mem_req", "work_dir", "scripts",  "accessKey", "access-key", "secretKey", "secret-key", "region",  "big_time",  'list_phenogc', 'cojo_slct_other', "paintor_fileannot", "paintor_listfileannot", "caviarbf_avalue", "gwas_cat", "genes_file", "genes_file_ftp",   "AMI", "instanceType", "instance-type", "bootStorageSize","maxInstances", "max-instances", "sharedStorageMount", "shared-storage-mount",'queue', "data", "pheno", "covariates", "cojo_slct_other", "file_gwas"]
allowed_params=allowed_params_input
allowed_params_cores=["plink_cpus_req", "gcta_cpus_req", "fm_cpus_req",'max_plink_cores']
allowed_params+=allowed_params_cores
allowed_params_intother=["max_forks", "n_pop","cojo_wind", "cojo_top_snps_chro"]
allowed_params+=allowed_params_intother
allowed_params_bolother=["cojo_slct", "cojo_actual_geno","cojo_top_snps"]
allowed_params+=allowed_params_bolother
allowed_params_float=["cut_maf", "threshold_p", "cojo_p"]
allowed_params+=allowed_params_float
allowed_params_memory=["gcta_mem_req" , "plink_mem_req", "other_mem_req", "boot-storage-size"]
allowed_params+=allowed_params_memory
params_filegwas=[ "head_beta", "head_se", "head_A1", "head_A2", "head_freq", "head_chr", "head_bp", "head_rs", "head_pval", "head_n"]
allowed_params+=params_filegwas





params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
        println "Check $parm  ************** is it a valid parameter -- are you using one rather than two - signs or vice-versa";
        System.exit(-1)
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
params.cojo_p=0
params.threshold_p=0
params.cojo_wind=10000
params.cut_maf=0.01
params.gcta_mem_req="15GB"
params.plink_mem_req="6GB"
params.gcta_cpus_req = 1
params.big_time='100h'
params.data=""
params.pheno=""
params.pos_cond      = 0
params.pos_ref = 0
params.chro_cond      = ""



gcta_mem_req=params.gcta_mem_req
gcta_cpus_req = params.gcta_cpus_req+2
plink_mem_req = params.plink_mem_req
max_plink_cores = params.max_plink_cores
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

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}

println "\nTesting data : ${params.data}\n"
println "Testing for gwas file : ${params.file_gwas}\n"

checkmultiparam(params,allowed_params_input, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_memory, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
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

plink_format=Channel.create()
plink_cojo=Channel.create()
plink_other=Channel.create()
raw_src_ch.separate (plink_format, plink_cojo, plink_other) { a -> [ a, a, a] }


check2 = Channel.create()
chrolist=chrolist.flatMap { list_str -> list_str.split() }.tap ( check2)

data_listind=Channel.fromPath(params.data)
process getListInd{
   input :
     file data from data_listind
   output:
     file(keepout) into (filekeepformat,filekeepcojo)
   script :
     varpheno=params.pheno!=""? " --pheno ${params.pheno}" : ""
   keepout="list_ind.keep"
   """
   gcta_extindplk.py --data ${data}  $varpheno --out $keepout
   """
}


gwas_format = Channel.fromPath(params.file_gwas)

process clean_plinik {
   cpus params.max_plink_cores
   memory params.plink_mem_req
   input :
     path(keepind) from filekeepformat
     tuple path(bed), path(bim), path(fam) from plink_format
   output :
      tuple path(keepind), path("${out}.bed"), path("${out}.bim"), path("${out}.fam") into clean_plink 
   script :
     listpos=params.pos_ref+","+params.pos_cond
     baseplk=bed.baseName
     out=params.output+"_gt"
     headkeep=params.data!="" ? " --keep $keepind " : ""
     """
     gcta_cleanplink.py --bfile $baseplk --chr $params.chro_cond $headkeep --out $out --threads ${params.max_plink_cores} --pos_list $listpos
     """
}

process doFormatData{
   cpus params.max_plink_cores
   memory params.plink_mem_req
   input :
     tuple path(keepind), path(bed), path(bim), path(fam) from clean_plink 
      path(gwas) from gwas_format
   output :
        tuple path(keepind), path(bed), path(bim), path(fam), path(out) into gwas_chro_cojo,gwas_chro_topsnp
   script :
      out=chro+".format"
      baseplk=bed.baseName
      headfreq=params.head_freq!="" ? " --freq_header ${params.head_freq}" : ""
      headn=params.head_n!="" ? " --freq_header ${params.head_n}" : ""
      headkeep=params.data!="" ? " --keep $keepind " : ""
      listpos=params.pos_ref+","+params.pos_cond
      """
      gcta_format.py --inp_asso $gwas  --rs_header ${params.head_rs} --pval_header ${params.head_pval} $headfreq --a1_header ${params.head_A1} --a2_header ${params.head_A2} --se_header ${params.head_se} --beta_header ${params.head_beta} --chro_header ${params.head_chr} --chr $chro --bfile $baseplk --out $out --threads ${params.max_plink_cores} --pos_list $listpos --updaters 1
      """
}


