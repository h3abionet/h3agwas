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

allowed_params = ['input_config', 'metal', 'gwama', 'heterogenity', 'metal_bin', 'GWAMA_bin', 'input_']

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

// external binary definition
params.metal_bin='metal'
params.gwama_bin='GWAMA'

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

info_file=configfile_analysis(params.file_config)
//liste_filesi_ch=Channel.fromPath(info_file[0]).merge(Channel.from(info_file[1]))
pos_file_ref=info_file[2]
if(pos_file_ref==-1){
println "not file reference foud see in config file column IsRefFile"
exit(0)
}
info_ref_rs=Channel.from(info_file[1][NumRef])
file_info_ref_rs=Channel.fromPath(info_file[0][NumRef])
process GetRsFile{
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
liste_filesi_ch=Channel.fromPath(info_file[0]).merge(Channel.from(info_file[1])).combine(file_rs_ref_chan)


process ChangeFormatFile {
    input :
      set file(file_assoc), val(info_file), file(file_ref) from liste_filesi_ch
    output :
      file(newfile_assoc) into (liste_file_gwama, liste_file_metal)
    script :
       newfile_assoc=file_assoc+".modif"
       """
       ma_change_format.py  --input_file $file_assoc --out_file $newfile_assoc --info_file $info_file --rs_ref $file_ref --sep_out TAB
       """
}

//(liste_file_gwama, liste_file_metal) = liste_files_mod_ch.toList().separate(2) { a -> [a, a] }

//liste_file_gwama = Channel.create()
//liste_file_metal = Channel.create()

//liste_files_mod_ch.separate(liste_file_gwama,liste_file_metal)

//liste_files_mod_ch=liste_files_mod_ch.collect()
liste_file_gwama=liste_file_gwama.collect()
liste_file_metal=liste_file_metal.collect()


if(params.gwama==1){
  //config channel
  process doGWAMA {
    input :
      val(list_file) from liste_file_gwama
    publishDir "${params.output_dir}/gwama", overwrite:true, mode:'copy'
    output :
      file("${out}*")
      //file("${out}.out") into res_gwama
    script :
      out = "gwama_res.stat"
      gc =  (params.ma_genomic_cont==1) ? "-gc -gco " : ""
      lfile=list_file.join(" ")
      """
      echo $lfile |awk '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' > fileListe
      ${params.gwama_bin}  --filelist fileListe --output $out -r $gc -qt --name_marker RSID --name_strand  DIRECTION --name_n N --name_eaf FREQA1 --name_beta BETA --name_se SE --name_ea A1 --name_nea A2
      """

  }
}

if(params.metal==1){

  process doMetal {
    input :
      val(list_file) from liste_file_metal
    publishDir "${params.output_dir}/metal", overwrite:true, mode:'copy'
    output :
      //file("${out}*") into res_metal
      file("${out}*") 
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

}
