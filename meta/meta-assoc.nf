#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Shaun Aron
 *      Rob Clucas
 *      Eugene de Beste
 *      Lerato Magosi
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS and metanalysis.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths

def helps = [ 'help' : 'help' ]

allowed_params = ['input_config', 'metal', 'gwama', 'heterogenity', 'metal_bin', 'GWAMA_bin', "ma_mem_req", "gwama_mem_req", "metasoft_mem_req", "plink_mem_req"]

/*
params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}
*/
def params_help = new LinkedHashMap(helps)


params.input_config=''
// what analysis
params.metal=0
params.gwama=0
params.metasoft=0
params.plink=0
params.mrmega=0
params.pheno="notdefined"

// external binary definition
params.metal_bin='metal'
params.gwama_bin='GWAMA'
params.mrmega_bin='MR-MEGA'
params.metasoft_bin="/opt/bin/Metasoft.jar"
params.plink_bin="plink"

params.ma_mrmega_pc=2
params.ma_random_effect=1
params.ma_genomic_cont=0
params.ma_inv_var_weigth=0
params.ma_mem_req="10G"
params.gwama_mem_req="20G"
params.metasoft_mem_req="20G"
params.plink_mem_req="10G"
params.big_time = '1000h'
params.file_config=''

params.metasoft_pvalue_table="/opt/bin/HanEskinPvalueTable.txt"
params.ma_metasoft_opt=""

params.ma_mrmega_opt=""

params.covariates=""
report_ch = Channel.empty()

gwama_mem_req=params.gwama_mem_req
ma_mem_req=params.ma_mem_req
metasoft_mem_req=params.metasoft_mem_req

def configfile_analysis(file){
   sep=','
   File theInfoFile = new File( file )
   if( !theInfoFile.exists() ) {
      println "File "+ file+" does not exist"
      exit(0)
   } 
   def lines = theInfoFile.readLines()
   def SplH=lines[0].split(sep)
   resInfo=[]
   resFile=[]
   lines.remove(0)
   PosCmp=-1
   CmtL=0
   for(line in lines){
       def splLine=line.split(sep)
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
       CmtL+=1
   }
 return([resFile,resInfo, NumRef])
}
checkexi=Channel.fromPath(params.file_config,checkIfExists:true)
info_file=configfile_analysis(params.file_config)
pos_file_ref=info_file[2]
if(pos_file_ref==-1){
println "not file reference foud see in config file column IsRefFile"
exit(0)
}
info_ref_rs=Channel.from(info_file[1][NumRef])
file_info_ref_rs=Channel.fromPath(info_file[0][NumRef],checkIfExists:true)
process GetRsFile{
    memory ma_mem_req
    input :
       file(file_assoc_rs) from file_info_ref_rs
       val(info_rs) from info_ref_rs
    output  :
       file(file_rs) into file_rs_ref_chan
    script :
        file_rs = file_assoc_rs+".rs"
        """
        ma_extract_rsid.py --input_file $file_assoc_rs --out_file $file_rs --info_file $info_rs 
        """
}
liste_filesi_ch=Channel.fromPath(info_file[0],checkIfExists:true).merge(Channel.from(info_file[1])).combine(file_rs_ref_chan)


process ChangeFormatFile {
    memory ma_mem_req
    input :
      set file(file_assoc), val(info_file), file(file_ref) from liste_filesi_ch
    output :
      file(newfile_assoc) into (liste_file_gwama, liste_file_metal, liste_file_metasoft, liste_file_mrmega)
      file(newfile_assocplk) into liste_file_plk
    script :
       newfile_assoc=file_assoc+".modif"
       newfile_assocplk=file_assoc+".modif.plk"
       """
       ma_change_format.py  --input_file $file_assoc --out_file $newfile_assoc --info_file $info_file --rs_ref $file_ref --sep_out TAB
       """
}

liste_file_gwama=liste_file_gwama.collect()
liste_file_metal=liste_file_metal.collect()
liste_file_metasoft=liste_file_metasoft.collect()
liste_file_mrmega=liste_file_mrmega.collect()
liste_file_plk=liste_file_plk.collect()


if(params.gwama==1){
 //config channel
  process doGWAMA {
     label 'metaanalyse'
     memory gwama_mem_req
     time params.big_time
    input :
      val(list_file) from liste_file_gwama
    publishDir "${params.output_dir}/gwama", overwrite:true, mode:'copy'
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
    publishDir params.output_dir, overwrite:true, mode:'copy'
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
    memory ma_mem_req
    time params.big_time
    input :
      val(list_file) from liste_file_mrmega
    publishDir "${params.output_dir}/mrmega", overwrite:true, mode:'copy'
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
      #echo $lfile |awk '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' > fileListe
      ma_printlistbonfile.py fileListe $lfile
      ${params.mrmega_bin}  --filelist fileListe -o $out $gc --qt --name_marker RSID --name_chr CHRO --name_pos POS --name_strand  DIRECTION --name_n N --name_eaf FREQA1 --name_beta BETA --name_se SE --name_ea A1 --name_nea A2 --pc ${params.ma_mrmega_pc} ${params.ma_mrmega_opt} --no_std_names
      """
  }
  process showMRMEGA {
    memory ma_mem_req
    publishDir params.output_dir
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
      val(list_file) from liste_file_metal
    publishDir "${params.output_dir}/metal", overwrite:true, mode:'copy'
    output :
      file("${out}1.stat") into res_metal
      set file("${out}1.stat"), file("${out}1.stat.info")
    script :
      out = "metal_res"
      lfile=list_file.join("\t")
      metal_config="metal_config.config"
      gc =  (params.ma_genomic_cont==1) ? " --genomic_control T " : "--genomic_control F"
      vw =  (params.ma_inv_var_weigth==1) ? " --inv_var_weigth T " : "--inv_var_weigth F"
      """
      echo $lfile |awk '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' > fileListe
      ma_get_configmetal.py --filelist fileListe  --output_configmetal $metal_config  $gc  $vw --out_file_metal $out
      ${params.metal_bin} $metal_config
      """
  }
  process showMetal {
    memory ma_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      file(assoc) from res_metal
    output:
      file("${out}*")  into report_Metal
    script:
      out = "metal"
      """
      metaanalyse_man.py  --inp $assoc --out ${out} --rs_header MarkerName --pval_header "P-value" --beta_header "Effect" --info_prog "Metal"
      """
  }
  report_ch = report_ch.flatten().mix(report_Metal.flatten())


}

if(params.metasoft==1){

  file_pvaltab=params.metasoft_pvalue_table
  process doMetaSoft {
    time params.big_time
    memory metasoft_mem_req
    label 'metaanalyse'
    input :
      val(list_file) from liste_file_metasoft
      //file(file_pvaltab) from filepvaltable
    publishDir "${params.output_dir}/metasoft", overwrite:true, mode:'copy'
    output :
      set file("${out}.meta"),file("${out}.res"),file("${out}.log"), file("${out}.files"), file("${out}.format.res"), file("${out}.pivot")
      file("${out}.format.res")  into res_metasoft
    script :
      out = "metasoft_res"
      lfile=list_file.join(" ")
      """
      ma_formatmetasoft.py $out $lfile
      java -jar ${params.metasoft_bin} -input $out".meta"  -output $out".res"   -log $out".log" -pvalue_table $file_pvaltab ${params.ma_metasoft_opt}
      ma_trans_outsetasoft.py $out".res" $out".files"  $out".format.res"
      """
  }
  process showMetasoft {
    time params.big_time
    memory metasoft_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      file(assoc) from res_metasoft
    output:
      file("${out}*")  into report_Metasoft
    script:
      out = "metasoft"
      """
      metaanalyse_man.py  --inp $assoc --out ${out} --rs_header RSID --pval_header "PVALUE_RE" --beta_header "BETA_RE" --info_prog "MetaSoft (Han and Eskin Random Effects model)"
      """
  }
  report_ch = report_ch.flatten().mix(report_Metasoft.flatten())

}

if(params.plink==1){
  process doPlinkMeta{
     time params.big_time
     memory params.plink_mem_req
     input :
      file(listeplk) from liste_file_plk
    publishDir "${params.output_dir}/plink", overwrite:true, mode:'copy'
     output :
       file("$out*") 
       file("${out}.meta") into  res_plink
    script :
     lpk=listeplk.join(" ")
     out=params.output+'_plink'
     """
     ${params.plink_bin} --meta-analysis $lpk + qt -out $out
     """

  }

  process showPlink {
    time params.big_time
    memory metasoft_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
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
  publishDir params.output_dir, overwrite:true, mode:'copy'
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




