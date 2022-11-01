#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2022
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS and metanalysis.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths

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



def helps = [ 'help' : 'help' ]

//allowed_params = ['input_config', 'metal', 'gwama', 'heterogenity', 'metal_bin', 'GWAMA_bin', "ma_mem_req", "gwama_mem_req", "metasoft_mem_req", "plink_mem_req"]
allowed_params_input = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts",  "high_ld_regions_fname", "accessKey", "access-key", "secretKey", "secret-key", "region",  "AMI", "instance-type", "boot-storage-size", "sharedStorageMount", "instanceType", "input_config","ma_metasoft_opt","metasoft_pvalue_table", "big_time", "file_config", "bootStorageSize","maxInstances","max-instances", "pheno","shared-storage-mount","csv_sep", "ma_mrmega_opt"]
allowed_params=allowed_params_input
allowed_params_cores=["max_plink_cores"]
allowed_params+=allowed_params_cores
allowed_params_intother=["ma_mrmega_pc"]
allowed_params+=allowed_params_intother
allowed_params_bin=["GWAMA_bin", "plink_bin", "metal_bin", "gwama_bin", "mrmega_bin", "metasoft_bin"]
allowed_params+=allowed_params_bin
allowed_params_bolother=["ma_random_effect", "ma_genomic_cont", "ma_inv_var_weigth","use_rs","ma_overlap_sample"]
allowed_params+=allowed_params_bolother
allowed_params_float=[]
allowed_params+=allowed_params_float
allowed_params_memory=["mrmega_mem_req", "ma_mem_req", "gwama_mem_req","metasoft_mem_req","plink_mem_req","ma_mem_req_utils"]
allowed_params+=allowed_params_memory
allowed_params_test=["plink", "metal", "gwama", "mrmega", "metasoft"]
allowed_params+=allowed_params_test



def params_help = new LinkedHashMap(helps)


params.input_config=''
// what analysis
params.metal=0
params.gwama=0
params.metasoft=0
params.csv_sep=','
params.plink=0
params.mrmega=0
params.pheno="meta"

// external binary definition
params.metal_bin='metal'
params.gwama_bin='GWAMA'
params.mrmega_bin='MR-MEGA'
params.metasoft_bin="/opt/bin/Metasoft.jar"
params.plink_bin="plink"
params.max_plink_cores=4


params.ma_mrmega_pc=2
params.ma_random_effect=1
params.mrmega_mem_req="20GB"
params.ma_genomic_cont=0
params.ma_inv_var_weigth=0
params.ma_mem_req="10G"
params.gwama_mem_req="20G"
params.metasoft_mem_req="20G"
params.plink_mem_req="10G"
params.big_time = '1000h'
params.file_config=''
params.use_rs=0
params.ma_overlap_sample=0

params.output_dir="output/"

params.metasoft_pvalue_table="/opt/bin/HanEskinPvalueTable.txt"
params.ma_metasoft_opt=""

params.ma_mrmega_opt=""

params.covariates=""
params.ma_mem_req_utils="10GB"
report_ch = Channel.empty()

gwama_mem_req=params.gwama_mem_req
ma_mem_req=params.ma_mem_req
metasoft_mem_req=params.metasoft_mem_req

params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}
def configfile_analysis(file, sep){
   File theInfoFile = new File( file )
   if(file.contains("s3://")){
     println "File " + file + "is in S3 so its existence can't be checked."
   }
   if(file.contains("az://")){
     println "File " + file + "is in Azure so its existence can't be checked."
   }
   if( !theInfoFile.exists() && !file.contains("s3://") && !file.contains("az://") ) {
      println "File "+ file+" does not exist"
      exit(0)
   } 
   def lines = theInfoFile.readLines()
   def SplH=lines[0].split(sep)
   if(SplH.size()==1){
      println SplH
      println "problem of tiy separator"
      exit(1)
   }
   resInfo=[]
   resFile=[]
   lines.remove(0)
   PosCmp=-1
   CmtL=0
   listnum=[]
   NumRef=-1
   for(line in lines){
      def splLine=line.split(sep)
      if(splLine.size()>4){
       def cmtelem=0
       def SubRes=[]
       while (cmtelem < SplH.size()){
            if(SplH[cmtelem]=='File'){
               resFile.add(splLine[cmtelem])
            }
            else if(SplH[cmtelem]=='IsRefFile' && splLine[cmtelem]=='1'){
               NumRef=CmtL
               cmtelem2=0
            }
            if(splLine[cmtelem]!='NA' && splLine[cmtelem]!='' && SplH[cmtelem]!='File' && SplH[cmtelem]!='IsRefFile'){
                 SubRes.add(SplH[cmtelem]+':'+splLine[cmtelem])
            }
            cmtelem+=1
       }       
       resInfo.add(SubRes.join(','))
       listnum.add(CmtL)
       CmtL+=1
     }
   }
 if(CmtL<2){
   println "no enough file found for meta analyse, check your column header and sep in your csv file\n exit"
   exit(1)
 }
 return([resFile,resInfo, NumRef, listnum])
}


checkmultiparam(params,allowed_params_input, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_memory, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_cores, java.lang.Integer, min=1, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_intother, java.lang.Integer, min=0, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_bolother, java.lang.Integer, min=0, max=null, possibleval=[0,1], notpossibleval=null)
checkmultiparam(params,allowed_params_test, java.lang.Integer, min=0, max=null, possibleval=[0,1], notpossibleval=null)
checkmultiparam(params,allowed_params_float, [java.lang.Float, java.lang.Integer, java.math.BigDecimal,java.lang.Double], min=0, max=null, possibleval=null, notpossibleval=null)


if(params.file_config==''){
println "not file config defined\n exit"
System.exit(-1)
}
checkexi=Channel.fromPath(params.file_config,checkIfExists:true)
info_file=configfile_analysis(params.file_config,params.csv_sep)
pos_file_ref=info_file[2]
liste_filesforref_ch=Channel.fromPath(info_file[0],checkIfExists:true).merge(Channel.from(info_file[1])).merge(Channel.from(info_file[3]))
process extract_rs_file{
    memory ma_mem_req
    input :
      set file(file_assoc), val(info_file), val(num) from liste_filesforref_ch
    publishDir "${params.output_dir}/rsinfo/", mode:'copy'
    output  :
       //file(file_rs) into file_rs_ref_chan, file_ref_rs_metal
       file(file_rs) into file_rs_formerge
    script :
        file_rs = num+"_"+file_assoc+".rs"
        """
        ma_extract_rsid.py --input_file $file_assoc --out_file $file_rs --info_file $info_file
        """
}

file_rs_formerge_m=file_rs_formerge.collect()

process extract_allrs{
    memory params.ma_mem_req_utils
   input :
     path(filers)  from file_rs_formerge_m
   publishDir "${params.output_dir}/rsinfo/", mode:'copy'
   output :
     path(listrs) into file_rs_ref_chan, file_ref_rs_metal
   script :
      listrs="all_rs"
      listfilers=filers.join(",")
      """  
      merge_rs.py --list_file $listfilers --out $listrs --use_rs ${params.use_rs}
      """
}
liste_filesi_ch=Channel.fromPath(info_file[0],checkIfExists:true).merge(Channel.from(info_file[1])).merge(Channel.from(info_file[3])).combine(file_rs_ref_chan)


/*deletedd file_ref*/
process ChangeFormatFile {
    memory ma_mem_req
    input :
      set file(file_assoc), val(info_file), val(num),file(file_ref) from liste_filesi_ch
    output :
      file(newfile_assoc) into (liste_file_gwama, liste_file_metal, liste_file_metasoft, liste_file_mrmega, liste_file_metasoft_extractinfo)
      file(newfile_assocplk) into liste_file_plk
    script :
       newfile_assoc="${file_assoc}_${num}.modif"
       newfile_assocplk="${file_assoc}_${num}.modif.plk"
       """
       ma_change_format_v2.py  --input_file $file_assoc --out_file $newfile_assoc --info_file $info_file  --sep_out TAB --use_rs ${params.use_rs} --rs_ref $file_ref
       """
}

liste_file_gwama=liste_file_gwama.collect()
liste_file_metal=liste_file_metal.collect()
liste_file_metasoft=liste_file_metasoft.collect()
liste_file_metasoft_extractinfo=liste_file_metasoft_extractinfo.collect()
liste_file_mrmega=liste_file_mrmega.collect()
liste_file_plk=liste_file_plk.collect()


if(params.gwama==1){
 //config channel
  process doGWAMA {
     label 'metaanalyse'
     memory gwama_mem_req
     time params.big_time
    input :
      file(list_file) from liste_file_gwama
    publishDir "${params.output_dir}/gwama", mode:'copy'
    output :
      //file("${out}*")
      file("${out}.out") into res_gwama
      set file("${out}.err.out"), file("${out}.gc.out"), file("${out}.log.out"), file("${out}.out")
    script :
      out = "gwama_res"
      gc =  (params.ma_genomic_cont==1) ? "-gc -gco " : ""
      optrandom = (params.ma_random_effect==1) ? " -r " : ""
      lfile=list_file.join(" ")
      """
      echo $lfile |awk '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' > fileListe
      ${params.gwama_bin}  --filelist fileListe --output $out $optrandom $gc -qt --name_marker RSID --name_strand  DIRECTION --name_n N --name_eaf FREQA1 --name_beta BETA --name_se SE --name_ea A1 --name_nea A2
      """
  }
  process showGWAMA {
    memory ma_mem_req
    publishDir params.output_dir, mode:'copy'
    input:
      file(assoc) from res_gwama
    output:
      file("${out}*")  into report_gwama
    script:
      out = "gwama"
      """
      metaanalyse_man.py  --inp $assoc --out ${out} --rs_header "rs_number" --pval_header "p-value" --beta_header "beta" --info_prog "GWAMA"
      """
  }
  report_ch = report_ch.flatten().mix(report_gwama.flatten())

}

if(params.mrmega==1){
  //config channel
  process doMRMEGA {
    label 'metaanalyse'
    memory params.mrmega_mem_req
    time params.big_time
    input :
      file(list_file) from liste_file_mrmega
    publishDir "${params.output_dir}/mrmega", mode:'copy'
    output :
      //file("${out}*")
      file("${out}.result") into res_mrmega
      //mrmega_res.log  mrmega_res.result
      set file("${out}.result"), file("${out}.log")
    script :
      out = "mrmega_res"
      gc =  (params.ma_genomic_cont==1) ? "--gc " : ""
      lfile=list_file.join(" ")
      """
      ma_printlistbonfile.py fileListe $lfile
      ${params.mrmega_bin}  --filelist fileListe -o $out $gc --qt --name_marker RSID --name_chr CHRO --name_pos POS --name_strand  DIRECTION --name_n N --name_eaf FREQA1 --name_beta BETA --name_se SE --name_ea A1 --name_nea A2 --pc ${params.ma_mrmega_pc} ${params.ma_mrmega_opt} --no_std_names
      """
  }
  process showMRMEGA {
    memory ma_mem_req
    publishDir "${params.output_dir}/mrmega", mode:'copy'
    input:
      file(assoc) from res_mrmega
    output:
      file("${out}*")  into report_mrmega
    script:
      out = "mrmega"
      """
      metaanalyse_man.py  --inp $assoc --out ${out} --rs_header "MarkerName" --pval_header "P-value_association" --beta_header "lnBF" --info_prog "MRMEGA"
      """
  }
  report_ch = report_ch.flatten().mix(report_mrmega.flatten())

}



if(params.metal==1){

  process doMetal {
    time params.big_time
    label 'metaanalyse'
    memory ma_mem_req
    input :
      file(list_file) from liste_file_metal
      file(file_ref_rs) from file_ref_rs_metal
      
    publishDir "${params.output_dir}/metal",  mode:'copy'
    output :
      file("${out}1.stat") into res_metal
      tuple file("${out}1.stat"), file("${out}1.stat.info")
      file("${out}_metal.format")
      file("${out}.log")
      file("${metal_config}")
    script :
      out = "metal_res"
      lfile=list_file.join("\t")
      metal_config="metal_config.config"
      gc =  (params.ma_genomic_cont==1) ? " --genomic_control T " : "--genomic_control F"
      vw =  (params.ma_inv_var_weigth==1) ? " --inv_var_weigth T " : "--inv_var_weigth F"
      sov=(params.ma_overlap_sample==1) ? " --overlap T " : " --overlap F "
      """
      echo $lfile |awk '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' > fileListe
      ma_get_configmetal.py --filelist fileListe  --output_configmetal $metal_config  $gc  $vw --out_file_metal $out $sov
      ${params.metal_bin} $metal_config &> ${out}.log
      merge_summarystat.py --input_file ${out}1.stat --info_file $file_ref_rs --out_file $out"_metal.format"
      """
  }
  process showMetal {
    memory ma_mem_req
    publishDir "${params.output_dir}/metal",  mode:'copy'
    input:
      file(assoc) from res_metal
    output:
      file("${out}*")  into report_Metal
    script:
      out = "metal"
      """
      metaanalyse_man.py  --inp $assoc --out ${out} --rs_header MarkerName --pval_header "P-value" --beta_header Zscore,Effect --info_prog "Metal"
      """
  }
  report_ch = report_ch.flatten().mix(report_Metal.flatten())


}

if(params.metasoft==1){

  //filepvaltable=channel.fromPath(params.metasoft_pvalue_table, checkIfExists:true)
  file_pvaltab=params.metasoft_pvalue_table
  process formatMetasoft{
    time params.big_time
    memory metasoft_mem_req
    input :
      file(list_file) from liste_file_metasoft
    publishDir "${params.output_dir}/metasoft",mode:'copy'
    output :
      set file("${out}.meta"), file("${out}.files"), file("${out}.pivot") into (format_sumstat, format_sumstat2)
    script :
      out = "metasoft_res"
      lfile=list_file.join(" ")
      """
      ma_formatmetasoft.py $out $lfile
      """
  }
  process extractInfoMetasoft{
    time params.big_time
    memory metasoft_mem_req
    input :
      file(list_file) from liste_file_metasoft_extractinfo
    publishDir "${params.output_dir}/metasoft",mode:'copy'
    output :
      file("${out}.N") into infoMetasoft
    script :
      out = "metasoft_res"
      lfile=list_file.join(" ")
      """
      ma_formatmetasoft_info.py $out $lfile
      """
  }
  process doMetasoft {
    time params.big_time
    memory metasoft_mem_req
    label 'metaanalyse'
    input :
      set file(metastat), file(listfiles), file(pivot) from format_sumstat
    publishDir "${params.output_dir}/metasoft",mode:'copy'
    output :
      set file("${out}.res"), file("${out}.log") into resbrut_metasoft
    script :
      out = "metasoft_res"
      """
      java -jar ${params.metasoft_bin} -input $metastat  -output $out".res"   -log $out".log" -pvalue_table $file_pvaltab ${params.ma_metasoft_opt}
      """
  }
  process FormatOutputMetasoft{
    time params.big_time
    memory metasoft_mem_req
    input :
      set path(metastat), path(listfiles), path(pivot) from format_sumstat2
      set path(metasoftres), path(filelog) from resbrut_metasoft
      path(InfoN) from infoMetasoft
    publishDir "${params.output_dir}/metasoft",mode:'copy'
    output :
      path("${out}.format.res") into res_metasoft
    script :
      out = "metasoft_res"
      """
      ma_trans_outsetasoft.py $metasoftres $listfiles $pivot  $InfoN $out".format.res"
      """
  }
  process showMetasoft {
    label 'metaanalyse'
    time params.big_time
    memory metasoft_mem_req
    publishDir "${params.output_dir}/metasoft",mode:'copy'
    input:
      file(assoc) from res_metasoft
    output:
      file("${out}*")  into report_Metasoft
    script:
      out = "metasoft"
      """
      metaanalyse_man.py  --inp $assoc --out ${out} --rs_header rs --pval_header "PVALUE_RE" --beta_header "BETA_RE" --info_prog "MetaSoft (Han and Eskin Random Effects model)"
      """
  }
  report_ch = report_ch.flatten().mix(report_Metasoft.flatten())

}

if(params.plink==1){
  process doPlinkMeta{
     time params.big_time
     memory params.plink_mem_req
     cpus params.max_plink_cores
     input :
      file(listeplk) from liste_file_plk
    publishDir "${params.output_dir}/plink", mode:'copy'
     output :
       file("$out*") 
       file("${out}.meta") into  res_plink
    script :
     lpk=listeplk.join(" ")
     out=params.output+'_plink'
     //weighted-z
     vw =  (params.ma_inv_var_weigth==1) ? " + weighted-z " : ""
     """
     ${params.plink_bin} --meta-analysis $lpk + qt $vw -out $out --threads ${params.max_plink_cores} 
     """

  }

  process showPlink {
    time params.big_time
    memory ma_mem_req
    publishDir "${params.output_dir}/plink", mode:'copy'
    input:
      file(assoc) from res_plink
    output:
      file("${out}*")  into report_plink
    script:
      out = "plink"
      """
      metaanalyse_man.py  --inp $assoc --out ${out} --rs_header SNP --pval_header "P" --beta_header "BETA" --info_prog "BETA"
      """
  }

   report_ch = report_ch.flatten().mix(report_plink.flatten())
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

nextflowversion =nextflow.version


if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"


//---- Modification of variables for pipeline -------------------------------//


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


process doReport {
  label 'latex'
  input:
    file(reports) from report_ch.toList()
  publishDir params.output_dir, mode:'copy'
  output:
    file("${out}.pdf")
  script:
    out = params.output+"-report"
    these_phenos     = params.pheno
    these_covariates = params.covariates
    config = getConfig()
    images = workflow.container
    texf   = "${out}.tex"
    template "make_assoc_report.py"
}


