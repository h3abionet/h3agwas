#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      jean-tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2019
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


def helps = [ 'help' : 'help' ]
allowed_params = ['file_listvcf', 'min_scoreinfo', "output_pat", "output_dir", "do_stat"]


params.help = false

nextflowversion =nextflow.version


if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"

report = new LinkedHashMap()

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}

params.queue    = 'batch'
params.remove_on_bp  = 1

params.file_listvcf=""
params.min_scoreinfo=0.6
params.max_plink_cores = 8 
//params.plink_mem_req = '10GB' // how much plink needs for this
params.other_mem_req = '10GB' // how much plink needs for this
params.output_pat="out"
params.output_dir="plink/"
params.statfreq_vcf="%AN %AC"
params.score_imp="INFO"
params.genetic_maps=""
params.do_stat=true
params.unzip_zip=0
params.unzip_password=""



if(params.file_listvcf==""){
error('params.file_listvcf : file contains list vcf not found')
}
/*read file to have list of channel for each vcf*/

if(params.unzip_zip){
list_vcf_unz=Channel.fromPath(file(params.file_listvcf).readLines())
process UnzipFile{
  input :
     file(vcfzip) from list_vcf_unz  
  publishDir "${params.output_dir}/vcfi", overwrite:true, mode:'copy'
  output :
   file("*.vcf.gz") into (list_vcf, list_vcf2)
   file("*.gz") 
  script :
    passd = (params.unzip_password!="") ? " -P ${params.unzip_password} " : ""
    """
    unzip $passd $vcfzip 
    """  
}

}else{
list_vcf=Channel.fromPath(file(params.file_listvcf).readLines())
list_vcf2=Channel.fromPath(file(params.file_listvcf).readLines())
}

if(params.do_stat){
process computedstat{
 memory params.plink_mem_req
  time   params.big_time
  input :
     file(vcf) from list_vcf2
  output :
     file("${Ent}") into listchrostat
  script :
    Ent=vcf.baseName+".stat"
    """
    bcftools query -f '%CHROM %REF %ALT %POS %INFO/${params.score_imp} ${params.statfreq_vcf}\n' $vcf > $Ent
    """
}
statmerg=listchrostat.collect()

process dostat{
 memory params.plink_mem_req
 input :
    file(allstat) from statmerg
 publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
 output :
   file("${fileout}*")
 script :
  fileout=params.output_pat+"_report"
  allfile=allstat.join(',') 
  """
  stat_vcf.py  --out $fileout --min_score ${params.min_scoreinfo} --list_files $allfile
  """
}
}
ref_ch=Channel.fromPath(params.reffasta)

if(params.min_scoreinfo>0){
list_vcf=ref_ch.combine(list_vcf)
process formatvcfscore{
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     set file(ref),file(vcf) from list_vcf
  output :
     set file("${Ent}.bed"), file("${Ent}.bim"), file("${Ent}.fam") into listchroplink
     file("${Ent}.bim") into listbimplink
  script :
    Ent=vcf.baseName
    """
    bcftools view -Ou -i '${params.score_imp}>${params.min_scoreinfo}' $vcf | bcftools norm -Ou -m -any | bcftools norm -Ou -f $ref |bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' |
  plink --bcf /dev/stdin \
    --keep-allele-order \
    --vcf-idspace-to _ \
    --const-fid \
    --allow-extra-chr 0 \
    --make-bed \
    --out ${Ent}
    cp ${Ent}.bim ${Ent}.save.bim
    awk \'{if(\$2==\".\"){\$2=\$1\":\"\$4\"_\"\$5\"_\"\$6};\$2=substr(\$2,1,20);print \$0}\' ${Ent}.save.bim > ${Ent}.bim
    """
}


}else{
process formatvcf{
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     file(vcf) from list_vcf
  output :
     set file("${Ent}.bed"), file("${Ent}.bim"), file("${Ent}.fam") into listchroplink
     file("${Ent}.bim") into listbimplink
  script :
    Ent=vcf.baseName
    """

    bcftools view -Ou $vcf | bcftools norm -Ou -m -any | bcftools norm -Ou -f $ref |bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' |
  plink --bcf /dev/stdin \
    --keep-allele-order \
    --vcf-idspace-to _ \
    --const-fid \
    --allow-extra-chr 0 \
    --make-bed \
    --out ${Ent}


    cp ${Ent}.bim ${Ent}.save.bim
    awk \'{if(\$2==\".\"){\$2=\$1\":\"\$4\"_\"\$5\"_\"\$6};\$2=substr(\$2,1,20);print \$0}\' ${Ent}.save.bim > ${Ent}.bim
    """
}
}
bimmerg=listbimplink.collect()
process GetRsDup{
    memory params.other_mem_req
    input :
      file(bim) from bimmerg
    output :
      set file(outdel),file(out) into duplicat
    script :
      lbim=bim.join(",")
      out="snpfile_red.rs"
      outdel="snpfile_del.rs"
      """
      search_dup_bim.py $lbim $out $outdel
      """ 
}
listchroplinkinf=duplicat.combine(listchroplink)
process TransformRsDup{
    input :
     set file(delrange),file(rstochange), file(bed),file(bim),file(fam) from listchroplinkinf
    output :
       set file("${newheader}.bed"),file("${newheader}.bim"),file("${newheader}.fam") into listchroplinkrs
       val(newheader) into plinkhead
    script :
       header=bed.baseName 
       newheader=header+"_rsqc"
       """
       cp $bim ${bim}.save
       replacers_forbim.py ${bim}.save $rstochange $bim  
       plink -bfile $header --keep-allele-order --make-bed -out $newheader --exclude range $delrange
       """
}
///dataE/AWIGenGWAS/shared/imputed_data_plink/Build/genetic_map_hg19.txt
if(params.genetic_maps!=""){
GMMap=Channel.fromPath(params.genetic_maps)
listchroplinkrsmap=GMMap.combine(listchroplinkrs)
process AddedCM{
    cpus params.max_plink_cores
    memory params.plink_mem_req
    time   params.big_time
    input :
       set map,file(bedi),file(bimi),file(fami) from listchroplinkrsmap
    output :
       set file(bedf),file(bimf),file(famf) into listchroplinkrsf
       val(header) into plinkheadf
    script :
       headeri=bedi.baseName
       header=headeri+"_map"
       bedf=header+".bed"
       bimf=header+".bim"
       famf=header+".fam"
       cm_shap=header+".shape"
       """
       chro=`head $bimi|awk '{print \$1}'|uniq`
       sed '1d' $map|awk -v chro=\$chro '{if(chro==\$1)print \$2"\\t"\$3"\\t"\$4}' >> $cm_shap
       awk '{print \$2}' $headeri".bim" | sort | uniq -d > duplicated_snps.snplist
       plink --bfile $headeri --exclude duplicated_snps.snplist --make-bed --keep-allele-order --out $headeri"_tmp"
       plink --bfile $headeri"_tmp" --list-duplicate-vars ids-only suppress-first
       plink --bfile $headeri"_tmp" --keep-allele-order --cm-map $cm_shap \$chro   --threads ${params.max_plink_cores} --make-bed --out $header  --exclude plink.dupvar
       """

}

}else{
plinkheadf=plinkhead
listchroplinkrsf=listchroplinkrs
}
listplinkinf=listchroplinkrsf.collect()
headplinkf=plinkheadf.collect()

process MergePlink{
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
       file(lplk) from listplinkinf
       val(hplk) from headplinkf
  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
  output :
     set file("${params.output_pat}.bed"), file("${params.output_pat}.bim"),file("${params.output_pat}.fam") into plinkformatf
  script :
       hplkFirst=hplk[0]
       hplk.remove(0)
       hplk2=hplk.join(',')
       """
       echo $hplk2|awk -F\',\' \'{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt"\\n"}\' > fileplk
       plink --bfile $hplkFirst --keep-allele-order --threads ${params.max_plink_cores} --merge-list fileplk --make-bed --out ${params.output_pat}
       """
}




