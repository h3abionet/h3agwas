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
nextflow.enable.dsl = 1



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
params.reffasta=""
params.data=''
params.genotype_pat="TYPED"


if(params.file_listvcf==""){
error('params.file_listvcf : file contains list vcf not found')
}
/*read file to have list of channel for each vcf*/

if(params.unzip_zip){
list_vcf_unz=Channel.fromPath(file(params.file_listvcf).readLines(), checkIfExists:true)
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
list_vcf=Channel.fromPath(file(params.file_listvcf).readLines(), checkIfExists:true)
list_vcf2=Channel.fromPath(file(params.file_listvcf).readLines())
}

if(params.do_stat){
process computedstat{
 label 'py3utils'
 memory params.plink_mem_req
  time   params.big_time
  input :
     file(vcf) from list_vcf2
  output :
     file("${Ent}") into listchrostat
  script :
    Ent=vcf.baseName+".stat"
    """
    bcftools query -f '%CHROM %REF %ALT %POS %INFO/$params.genotype_pat %INFO/${params.score_imp} ${params.statfreq_vcf}\n' $vcf > $Ent
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
   tuple path("${fileout}.tex"),path('*report_freq.pdf'), path("*_report_score.pdf") into tex_stat
 script :
  fileout=params.output_pat+"_report"
  allfile=allstat.join(',') 
  """
  stat_vcf_v2.py  --out $fileout --min_score ${params.min_scoreinfo} --list_files $allfile
  """
}
process texfile{
  label 'latex'
  input :
    tuple path(tex), path(freqpdf), path(scorepdf) from tex_stat
  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
  output :
   path("${out}.pdf")
  script :
     out=tex.baseName
     """
    pdflatex $out >& /dev/null
    pdflatex $out
    """ 
}
}
//ref_ch=Channel.fromPath(params.reffasta, checkIfExists:true)


hgrefi=Channel.fromPath(params.reffasta, checkIfExists:true)                    
process checkfasta{                                                             
  cpus params.max_plink_cores                                                   
  label 'py3utils'                                                              
  errorStrategy { task.exitStatus == 1 ? 'retry' : 'terminate' }                
  maxRetries 1                                                                  
  input :                                                                       
     path(fasta) from hgrefi                                                    
  output :                                                                      
     tuple path("$fasta2"), path("${fasta2}.fai") into ref_ch
  script :                                                                      
    fasta2="newfasta.gz"                                                        
    """                                                                         
    if [ "${task.attempt}" -eq "1" ]                                            
    then                                                                        
    samtools faidx $fasta                                                       
    cp $fasta $fasta2                                                           
    mv $fasta".fai"  $fasta2".fai"                                              
    else                                                                        
    zcat $fasta | awk '{print \$1}' | bgzip -@ ${params.max_plink_cores} -c > $fasta2               
    samtools faidx $fasta2                                                      
    fi                                                                          
    """                                                                         
} 
list_vcf=ref_ch.combine(list_vcf)

if(params.min_scoreinfo>0){
process formatvcfscore{
  label 'py3utils'
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     set file(ref), file(refidx),file(vcf) from list_vcf
  publishDir "${params.output_dir}/bcf_filter", overwrite:true, mode:'copy', pattern: "*.bcf"
  output :
     set file("${Ent}.bed"), file("${Ent}.bim"), file("${Ent}.fam") into listchroplink
     file("${Ent}.bim") into listbimplink
     set val(Ent), file("${Ent}.bcf") into bcf_ch
  script :
    Ent=vcf.baseName.replaceAll(/.vcf$/, '')
    threadn=params.max_plink_cores/ 5
    threadn = threadn.round()
    """
    hostname
    bcftools view -Ou -i '${params.score_imp}>${params.min_scoreinfo}' $vcf --threads $threadn | bcftools norm -Ou -m -any --threads  $threadn | bcftools norm -Ou -f $ref --threads $threadn  --threads  $threadn |bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' --threads  $threadn > $Ent".bcf"
  cat $Ent".bcf" |plink --bcf /dev/stdin \
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

process convertbcf_invcf{
  label 'py3utils'
  time   params.big_time
  publishDir "${params.output_dir}/vcf_filter", overwrite:true, mode:'copy', pattern: "*.vcf.gz"
  input :
      set val(Ent), file(bcf) from bcf_ch
   output :
      val(vcf) 
   script :
     vcf=Ent+".vcf.gz"
     """
     bcftools convert $bcf -O z > $vcf
     """ 
}


}else{
process formatvcf{
  label 'py3utils'
  cpus params.max_plink_cores
  memory params.plink_mem_req
  time   params.big_time
  input :
     tuple path(ref),path(vcf) from list_vcf
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
if(params.genetic_maps!=""){
GMMap=Channel.fromPath(params.genetic_maps, checkIfExists:true)
listchroplinkrsmap=GMMap.combine(listchroplinkrs)
process AddedCM{
    cpus params.max_plink_cores
    memory params.plink_mem_req
    time   params.big_time
    input :
       set file(map),file(bedi),file(bimi),file(fami) from listchroplinkrsmap
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

if(params.data!=''){

data_ch=channel.fromPath(params.data, checkIfExists:true)
process update_name{
  label 'R'
  cpus params.max_plink_cores
  memory params.plink_mem_req
 input : 
   tuple path(bed), path(bim), path(fam) from plinkformatf
   path(data) from data_ch 
 output :
     tuple path("${params.output_pat}_idupdate.bed"), path("${params.output_pat}_idupdate.bim"),path("${params.output_pat}_idupdate.fam") into plinkformatf_chgid
     path("${out}.log")
 script :
     out=params.output_pat+"_idupdate"
     bfile=bed.baseName
     """
     change_updateidplink.r $fam $data update_id
     plink -bfile $bfile -make-bed --keep-allele-order --update-ids update_id -out $out
     """
}

}



