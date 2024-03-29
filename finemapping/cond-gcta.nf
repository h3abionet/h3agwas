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


allowed_params_input = ["input_dir","input_pat","output","output_dir","plink_mem_req", "work_dir", "scripts",  "accessKey", "access-key", "secretKey", "secret-key", "region",  "big_time",  'list_phenogc', 'cojo_slct_other', "paintor_fileannot", "paintor_listfileannot", "caviarbf_avalue", "gwas_cat", "genes_file", "genes_file_ftp",   "AMI", "instanceType", "instance-type", "bootStorageSize","maxInstances", "max-instances", "sharedStorageMount", "shared-storage-mount",'queue', "data", "pheno", "covariates", "cojo_slct_other", "file_gwas", 'gcta_bin']
allowed_params=allowed_params_input
allowed_params_cores=["plink_cpus_req", "gcta_cpus_req", "fm_cpus_req",'max_plink_cores']
allowed_params+=allowed_params_cores
allowed_params_intother=["max_forks", "n_pop", "pos_ref", 'around']
allowed_params+=allowed_params_intother
allowed_params_bolother=[ "cojo_actual_geno", "multi_cond"]
allowed_params+=allowed_params_bolother
allowed_params_float=["cut_maf", "threshold_p", "cojo_p"]
allowed_params+=allowed_params_float
allowed_params_memory=["gcta_mem_req" , "plink_mem_req", "other_mem_req", "boot-storage-size"]
allowed_params+=allowed_params_memory
params_filegwas=[ "head_beta", "head_se", "head_A1", "head_A2", "head_freq", "head_chr", "head_bp", "head_rs", "head_pval", "head_n"]
allowed_params+=params_filegwas
params_strorint=["chro_cond", 'pos_cond']
allowed_params+=params_strorint





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
params.cojo_actual_geno=0
params.around=100000
params.multi_cond=1



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
checkmultiparam(params,allowed_params_float, [java.lang.Double,java.lang.Float, java.lang.Integer, java.math.BigDecimal], min=0, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,params_filegwas, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,params_strorint, [java.lang.String,java.lang.Integer], min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_bolother, java.lang.Integer, min=0, max=null, possibleval=[0,1], notpossibleval=null)




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



data_listind=Channel.fromPath(params.data,checkIfExists:true)
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


gwas_format = Channel.fromPath(params.file_gwas,checkIfExists:true)

if(params.pos_ref<1){
 process clean_plinik {
   cpus params.max_plink_cores
   memory params.plink_mem_req
   input :
     path(keepind) from filekeepformat
     tuple path(bed), path(bim), path(fam) from plink_format
   output :
      tuple path(keepind), path("${out}.bed"), path("${out}.bim"), path("${out}.fam") into clean_plink, clean_plink2,plk_pos_ch
   script :
     listpos=params.pos_ref+","+params.pos_cond
     baseplk=bed.baseName
     out=params.output+"_gt"
     headkeep=params.data!="" ? " --keep $keepind " : ""
     chroC=params.chro_cond
     posCmin=params.pos_cond.toString().split(',').collect{it as int }.min()
     posCmax=params.pos_cond.toString().split(',').collect{it as int }.max()
     posCmin=(params.pos_ref==0) ? "$posCmin" : "${[posCmin, params.pos_ref].min()}"
     posCmax=(params.pos_ref==0) ? "$posCmax" : "${[posCmax, params.pos_ref].max()}"
     filerange="bedfile"
     """
     begin=`expr $posCmin - ${params.around}`
     end=`expr $posCmax + ${params.around}`
    if [ "\$begin" -lt 1 ]
     then
     begin=1
     fi
     echo -e "$chroC\t\$begin\t\$end\t$chroC:$posCmin" > $filerange
     plink --keep-allele-order --bfile $baseplk --make-bed --out $out --extract range $filerange --maf ${params.cut_maf} $headkeep
     cp   $out".bim" $out".save.bim"
     awk -F "\\t" 'OFS="\\t" {\$2=\$1":"\$4;print \$0}' $out".save.bim" > $out".bim"
     """
 }


 process format_data_wind{
   cpus params.max_plink_cores
   memory params.plink_mem_req
   input :
     tuple path(keepind), path(bed), path(bim), path(fam) from clean_plink 
      path(gwas) from gwas_format
   output :
        tuple path(keepind), path(bed), path(bim), path(fam), path(out) into (gwas_chro_cojo , gwas_chro_cojo_all)
   script :
      out=params.output+".format"
      baseplk=bed.baseName
      headfreq=params.head_freq!="" ? " --freq_header ${params.head_freq}" : ""
      headn=params.head_n!="" ? " --n_header ${params.head_n}" : ""
      headkeep=params.data!="" ? " --keep $keepind " : ""
     chroC=params.chro_cond
     posCmin=params.pos_cond.toString().split(',').collect{it as int }.min()
     posCmax=params.pos_cond.toString().split(',').collect{it as int }.max()
     posCmin=(params.pos_ref==0) ? "$posCmin" : "${[posCmin, params.pos_ref].min()}"
     posCmax=(params.pos_ref==0) ? "$posCmax" : "${[posCmax, params.pos_ref].max()}"
      """
      begin=`expr $posCmin - ${params.around}`
      end=`expr $posCmax + ${params.around}`
      gcta_format.py --inp_asso $gwas  --rs_header ${params.head_rs} --pval_header ${params.head_pval}  --a1_header ${params.head_A1} --a2_header ${params.head_A2} --se_header ${params.head_se} --beta_header ${params.head_beta} --chro_header ${params.head_chr} --chr $params.chro_cond --bfile $baseplk --out $out --threads ${params.max_plink_cores} --begin \$begin --end \$end --updaters 1 --bp_header ${params.head_bp} $headn $headfreq --print_pos 0
      """
 }
}else{
 process clean_plink_listpos {
   cpus params.max_plink_cores
   memory params.plink_mem_req
   input :
     path(keepind) from filekeepformat
     tuple path(bed), path(bim), path(fam) from plink_format
   output : 
      tuple path(keepind), path("${out}.bed"), path("${out}.bim"), path("${out}.fam") into clean_plink, clean_plink2,plk_pos_ch
   script :
     listpos=params.pos_ref+","+params.pos_cond
     baseplk=bed.baseName
     out=params.output+"_gt"
     headkeep=params.data!="" ? " --keep $keepind " : ""
     chroC=params.chro_cond
     listpos=params.pos_ref+','+params.pos_cond
     """
     gcta_cleanplink.py --out $out --chr $chroC --bfile $baseplk $headkeep --threads ${params.max_plink_cores} --pos_list $listpos
     """
 }
 process format_data_wind_listpos{
   cpus params.max_plink_cores
   memory params.plink_mem_req
   input :
     tuple path(keepind), path(bed), path(bim), path(fam) from clean_plink
      path(gwas) from gwas_format
   output :
        tuple path(keepind), path(bed), path(bim), path(fam), path(out) into gwas_chro_cojo,gwas_chro_cojo_all
   script :
      out=params.output+".format"
      baseplk=bed.baseName
      headfreq=params.head_freq!="" ? " --freq_header ${params.head_freq}" : ""
      headn=params.head_n!="" ? " --n_header ${params.head_n}" : ""
      headkeep=params.data!="" ? " --keep $keepind " : ""
      chroC=params.chro_cond
      listpos=params.pos_ref+','+params.pos_cond
      """
      gcta_format.py --inp_asso $gwas  --rs_header ${params.head_rs} --pval_header ${params.head_pval}  --a1_header ${params.head_A1} --a2_header ${params.head_A2} --se_header ${params.head_se} --beta_header ${params.head_beta} --chro_header ${params.head_chr} --chr $params.chro_cond --bfile $baseplk --out $out --threads ${params.max_plink_cores} --updaters 1 --bp_header ${params.head_bp} $headn $headfreq --print_pos 0 --list_pos $listpos
      """
 }

}
process computed_ld{
    input :
      tuple path(keepind), path(bed), path(bim), path(fam) from clean_plink2
    output :
      file(ld) into ld_res_notsq
      file(ld2) into ld_res_sq
    publishDir "${params.output_dir}/ld/",  mode:'copy'
    script :
       listpos= (params.pos_ref==0) ?"${params.pos_cond} "  : "${params.pos_cond},${params.pos_ref}"
       chr=params.chro_cond
       base        =  bed.baseName
       out = params.output
       ld=params.output+'_notsq.ld'
       ld2=params.output+'_sq.ld'
       """
       echo $listpos|awk -v chro=$chr -F"," '{for(Cmt=1;Cmt<=NF;Cmt++)print chro"\t"\$Cmt"\t"\$Cmt"\t"chro":"\$Cmt}' >ldlist
       plink -bfile $base --r2           --ld-window 99999  --out ${params.output}_notsq    --ld-window-kb 10000 --ld-window-r2 0
       plink -bfile $base --r2 square --out ${params.output}_sq
       """
}

/*
process plot_ld{
   label 'R'
  input :
     file(ld2) from ld_res_sq
     set path(keepind),file(bim), file(bed), file(fam) from plk_pos_ch
  output :
     file(fileout) into plot_ld
  publishDir "${params.output_dir}/res/fig/",  mode:'copy'
  script :
      fileout=params.output+'_ld.pdf'
      posref= (params.pos_ref<1) ?""  : " --pos_ref  ${params.pos_ref}"
      """
      cond_plotld.r --ld $ld2 --bim $bim --out $fileout  $posref
      """
}*/

list_cond=channel.of(params.pos_cond.toString().split(','))
process run_condgcta {
   label 'gcta'
   cpus params.gcta_cpus_req
   memory params.gcta_mem_req
   input :
      tuple path(keepind), path(bed), path(bim), path(fam), path(gwas) from gwas_chro_cojo
   each poscond from list_cond 
   publishDir "${params.output_dir}/condgcta/",  mode:'copy'
   output :
     file("${out}.cma.cojo") into cond_pos
     file("$out*")
   script :    
      baseplk=bed.baseName
      chro=params.chro_cond
      out=poscond+'_cond'
      cojoactu=params.cojo_actual_geno==1 ? " --cojo-actual-geno " : ""
      """
      echo "$chro:$poscond" > listcond
      ${params.gcta_bin} --bfile $baseplk --chr $chro --out $out --cojo-file ${gwas}   --thread-num ${params.gcta_cpus_req} --cojo-cond listcond $cojoactu  --cojo-collinear 0.99
      """
}
cond_pos=cond_pos.collect()
if(params.pos_cond.toString().split(',').size()>1 & params.multi_cond==1){
process run_condgcta_all {
   label 'gcta'
   cpus params.gcta_cpus_req
   memory params.gcta_mem_req
   input :
      tuple path(keepind), path(bed), path(bim), path(fam), path(gwas) from gwas_chro_cojo_all
   publishDir "${params.output_dir}/condgcta/",  mode:'copy'
   output :
     file("${out}.cma.cojo") into cond_pos_all
     file("$out*")
   script :
      baseplk=bed.baseName
      chro=params.chro_cond
      out='all_cond'
      cojoactu=params.cojo_actual_geno==1 ? " --cojo-actual-geno " : ""
      """
      echo "${params.pos_cond}" |awk -F"," '{for(col=1;col<=NF;col++)print ${params.chro_cond}":"\$col}' > listcond
      ${params.gcta_bin} --bfile $baseplk --chr $chro --out $out --cojo-file ${gwas}   --thread-num ${params.gcta_cpus_req} --cojo-cond listcond $cojoactu  --cojo-collinear 0.99
      """
}
cond_pos=cond_pos.combine(cond_pos_all).collect()
}
process merge_condld{
   label 'R'
   input :
    path(allfile) from cond_pos 
    path(ld) from ld_res_notsq
   publishDir "${params.output_dir}/condgcta/",  mode:'copy'
   output :
      path("$out*")
   script :
     out=params.output+".csv"
     allfile=allfile.join(",")
     posref=(params.pos_ref>0)? " --pos_ref ${params.pos_ref} " : ""
     """
     merge_gctacond.r --ld $ld --files_cond $allfile --out $out $posref
     """
}
