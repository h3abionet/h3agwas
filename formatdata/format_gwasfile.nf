#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      jean-tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2022
 *
 *
 * Description  : Nextflow pipeline to transform vcf file in plink and other format
 *
 *(C) University of the Witwatersrand, Johannesburg, 2016-2019 on behalf of the H3ABioNet Consortium
 *This is licensed under the MIT Licence. See the "LICENSE" file for details
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;
nextflow.enable.dsl = 1

/*definition*/
def errormess(message,exitn=0){
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
    errormess(messageerror,2)
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




filescript=file(workflow.scriptFile)
projectdir="${filescript.getParent()}"
dummy_dir="${projectdir}/../qc/input"


// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}

def helps = [ 'help' : 'help' ]
allowed_params_file = ['file_gwas', 'file_ref_gzip', "output_dir","output"]
allowed_header = ['head_pval', 'head_freq', 'head_bp', 'head_chr', 'head_rs', 'head_beta', 'head_se', 'head_A1', 'head_A2', 'sep']
allowed_headernew = ['headnew_pval', 'headnew_freq', 'headnew_bp', 'headnew_chr', 'headnew_rs', 'headnew_beta', 'headnew_se', 'headnew_A1', 'headnew_A2']
optional_param=["input_dir", "input_pat"]
//chro_ps 0 --bp_ps 1 --rs_ps 
allowed_posfilref=['poshead_chro_inforef', 'poshead_bp_inforef','poshead_rs_inforef']
allowed_params=allowed_params_file
allowed_params+=allowed_header
allowed_params+=allowed_headernew



params.mem_req = '2GB' // how much plink needs for this
params.file_gwas=""

params.head_pval = ""
params.head_freq = ""
params.head_bp = ""
params.head_chr = ""
params.head_rs = ""
params.head_beta=""
params.head_se=""
params.head_A1=""
params.head_A2=""
params.head_N=""
params.sep=""

params.headnew_pval = ""
params.headnew_freq = ""
params.headnew_bp = ""
params.headnew_chr = ""
params.headnew_rs = ""
params.headnew_beta=""
params.headnew_se=""
params.headnew_A1=""
params.headnew_A2=""
params.headnew_N=""

params.poshead_chro_inforef=0
params.poshead_bp_inforef=1
params.poshead_rs_inforef=2
params.poshead_a1_inforef=3
params.poshead_a2_inforef=4

if(params.file_gwas==""){
error('params.file_gwas: file contains gwas not found')
}


checkmultiparam(params,allowed_params_file, java.lang.String, min=null, max=null, possibleval=null, notpossibleval="")
checkmultiparam(params,allowed_header, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_headernew, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,optional_param, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_posfilref, [java.lang.String,java.lang.Integer], min=null, max=null, possibleval=null, notpossibleval=null)

headnew_freq=params.headnew_freq
headnew_pval=params.headnew_pval
headnew_beta=params.headnew_beta
headnew_A1=params.headnew_A1
headnew_A2=params.headnew_A2
headnew_N=params.headnew_N
headnew_se=params.headnew_N


if(headnew_pval=="")headnew_pval=params.head_pval 
if(headnew_freq=="")headnew_freq=params.head_freq 
if(headnew_beta=="")headnew_beta=params.head_beta 
if(headnew_se=="")headnew_se=params.head_se 
if(headnew_A1=="")headnew_A1=params.head_A1 
if(headnew_A2=="")headnew_A2=params.head_A2 
//if(params.headnew_N=="")params.headnew_N=params.head_N 

 if(params.input_dir!="" || params.input_pat!=''){
 print("used plink file")
 bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
 bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
 fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
 }else{
 bed="${dummy_dir}/00"
 bim="${dummy_dir}/01"
 fam="${dummy_dir}/02"
 }
 fileplk= Channel.create()
 Channel
     .from(file(bed),file(bim),file(fam))
     .buffer(size:3)
     .map { a -> [a[0], a[1], a[2]] }
     .set { fileplk }


gwas_chrolist = Channel.fromPath(params.file_gwas, checkIfExists:true)

if(params.head_chr!=""){
 headnew_bp=params.headnew_bp
 headnew_chr=params.headnew_chr
 if(headnew_bp=="")headnew_bp=params.head_bp 
 if(headnew_chr=="")headnew_chr=params.head_chr 
 headnew_rs="rs"
 if(params.head_rs!="")headnew_rs=params.head_rs
 if(params.headnew_rs!="") headnew_rs=params.head_rs 
 gwas_chrolist_ext = Channel.fromPath(params.file_gwas)
 process getListeChro{
         input :
           file(gwas_res) from gwas_chrolist
         output :
           file("filechro") into chrolist
         script:
          sep=(params.sep!="") ?  " --sep ${params.sep}" : ""
          """
          extractlistchro.py --input_file $gwas_res --chro_header ${params.head_chr} $sep > filechro
         """
 }
 
 


 chrolist2=Channel.create()
 chrolist.flatMap { list_str -> list_str.readLines()[0].split() }.set { chrolist2 }
 
 process ExtractChroGWAS{
     memory params.mem_req 
     input :
       file(gwas) from gwas_chrolist_ext
     each chro from  chrolist2
     output :
       set val(chro), file(gwas_out) into gwas_format_chro
     script :
       gwas_out=gwas.baseName+"_"+chro+".gwas"
       sep=(params.sep!="") ?  " --sep ${params.sep}" : ""
       infofile="Chro:${params.head_chr}:${headnew_chr},Pos:${params.head_bp}:${headnew_bp},A2:${params.head_A2}:${headnew_A2},A1:${params.head_A1}:${headnew_A1},af:${params.head_freq}:${headnew_freq},Beta:${params.head_beta}:${headnew_beta},Se:${params.head_se}:${headnew_se},Pval:${params.head_pval}:${headnew_pval},N:${params.head_N}:${params.headnew_N},SNP:${params.head_rs}:${params.headnew_rs}"
       """
       extractandformat_gwas.py --input_file $gwas --out_file ${gwas_out} --chr $chro --info_file $infofile
       """
 }
 if(params.file_ref_gzip==""){
 error('params.file_ref_gzip : file contains information for rs notnot found')
 }
 gwas_format_chro_rs=gwas_format_chro.combine(Channel.fromPath(params.file_ref_gzip, checkIfExists:true))
 if(params.head_rs==""){ 
   process ExtractRsIDChro{
     memory params.mem_req 
     input :
      set val(chro), file(gwas), file(rsinfo) from gwas_format_chro_rs
     output :
       set val(chro), file(gwas),file(outrs) into rsinfo_chro
     script :
       outrs="info_rs_"+chro+".rs"
     """
     zcat $rsinfo | extractrsid_bypos.py --file_chrbp $gwas --out_file $outrs --ref_file stdin --chr $chro --chro_ps ${params.poshead_chro_inforef} --bp_ps ${params.poshead_bp_inforef} --rs_ps ${params.poshead_rs_inforef} --a1_ps ${params.poshead_a1_inforef}  --a2_ps ${params.poshead_a2_inforef} 
     """ 
  }
 }else{
  rsinfo_chro=gwas_format_chro.combine(channel.fromPath('${dummy_dir}/03') )
 }
 
 
 rsinfo_chroplk=rsinfo_chro.combine(fileplk)
 process MergeRsGwasChro{
     memory params.mem_req 
     input :
       set val(chro), file(gwas),file(chrors), file(bed), file(bim), file(fam) from rsinfo_chroplk
     output :
       file(outmerge) into gwas_rsmerge
     script :
      outmerge="merge_"+chro+".gwas"
      bfileopt= (params.input_pat!="" || params.input_dir!="") ?  " --bfile "+bed.baseName : ""
      Nheadopt=(params.head_N!="") ? " --N_head ${params.head_N} " : ""
      Freqheadopt=(params.head_freq!="") ? " --freq_head ${params.head_freq} " : ""
 
      NheadNewopt=(params.headnew_N!="") ? " --Nnew_head ${params.headnew_N} " : ""
      FreqNewheadopt=(headnew_freq!="") ? " --freqnew_head ${headnew_freq} " : ""
      addrsopt=(params.head_rs=="") ? " --input_rs $chrors " : ""
      """
      mergeforrs.py --input_gwas $gwas $addrsopt   --out_file $outmerge --chro_head  ${headnew_chr} --bp_head  ${headnew_bp} --rs_head ${headnew_rs} --chro $chro $bfileopt  $Nheadopt $Freqheadopt $NheadNewopt $FreqNewheadopt  --a1_head ${headnew_A1} --a2_head  ${headnew_A2}
      """
 
 }

 gwas_rsmerge_all=gwas_rsmerge.collect()
 
 process MergeAll{
    memory params.mem_req 
    input :
       file(allfile) from gwas_rsmerge_all
    publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
    output :
       file(fileout) 
    script :
      file1=allfile[0]
      listefiles=allfile.join(" ")
      fileout=params.output
      """
      head -1 $file1 > $fileout
      ls $listefiles  | xargs -n 1 tail -n +2 >> $fileout
      """
 
 }
 
}else{

 headnew_chr="chr"
 if(params.head_chr!="")headnew_chr=params.head_chr
 if(params.headnew_chr!="") headnew_chr=params.headnew_chr

 headnew_rs="rs"
 if(params.head_rs!="")headnew_rs=params.head_rs
 if(params.headnew_rs!="") headnew_rs=params.headnew_rs

 headnew_bp="bp"
 if(params.head_bp!="")headnew_bp=params.head_bp
 if(params.headnew_bp!="") headnew_bp=params.headnew_bp

 gwas_ch= Channel.fromPath(params.file_gwas, checkIfExists:true)
 gwas_ch2= Channel.fromPath(params.file_gwas, checkIfExists:true)
 rsinfo_ch=Channel.fromPath(params.file_ref_gzip,checkIfExists:true)

 process format_sumstat{
     memory params.mem_req
     input :
       file(gwas) from gwas_ch2
     output :
       file(gwas_out) into gwas_format
     script :
       gwas_out=gwas.baseName+"_tmp.gwas"
       sep=(params.sep!="") ?  "--sep ${params.sep}" : ""
       infofile="A2:${params.head_A2}:${headnew_A2},A1:${params.head_A1}:${headnew_A1},af:${params.head_freq}:${headnew_freq},Beta:${params.head_beta}:${headnew_beta},Se:${params.head_se}:${headnew_se},Pval:${params.head_pval}:${headnew_pval},N:${params.head_N}:${params.headnew_N},SNP:${params.head_rs}:${headnew_rs}"
       """
       extractandformat_gwas.py --input_file $gwas --out_file ${gwas_out}  --info_file $infofile
       """
 }

 process extract_rsid{
  memory params.mem_req
  input :
     file(gwas) from gwas_ch
     file(rsinfo) from rsinfo_ch
   output :
      file(outrs) into rsout_ch
   script:
      sep=(params.sep!="") ?  " --sep ${params.sep}" : ""
      outrs="chrbp_rs.out"
      """
     zcat $rsinfo | extractchrpos_byrs.py --gwas $gwas --out_file $outrs --ref_file stdin  --chro_ps ${params.poshead_chro_inforef} --bp_ps ${params.poshead_bp_inforef} --rs_ps ${params.poshead_rs_inforef} --a1_ps ${params.poshead_a1_inforef}  --a2_ps ${params.poshead_a2_inforef} $sep --rs_head ${params.head_rs}  --a1_head ${params.head_A1} --a2_head ${params.head_A2}
      """
 }
 process add_rsid{
   input :
      file(gwas) from gwas_format
      file(outrs) from rsout_ch
   publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
   output :
       file("$gwasf*")
   script :
     gwasf=params.output
     println gwasf
     """
     add_chrbp_byrsid.py --gwas  $gwas --chrbp_info $outrs --rs_head ${params.head_rs} --chronew_head ${headnew_chr} --bpnew_head ${headnew_bp} --out $gwasf
     """
 }

}



 
