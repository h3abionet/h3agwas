#!/usr/bin/env nextflow

process check_names_plkconvert{                                                 
   label 'R'                                                                    
   input :                                                                      
       tuple path(bed), path(bim), path(fam)                                    
       path(data)                                                               
   publishDir "${params.output_dir}/format/plink/",  mode:'copy'                
   output :                                                                     
       tuple path("${newplk}.bed"), path("${newplk}.bim"), path("${newplk}.fam")
   script :                                                                     
     plk=bed.baseName                                                           
     newplk=plk+'_idupdate'                                                     
     """                                                                        
     cp $bed $newplk".bed"                                                      
     cp $bim $newplk".bim"                                                      
     change_names_plkconvert.r $fam $data $newplk".fam"                         
     """                                                                        
}      

process unzip_folder{
  input :
     path(vcfzip), path(outputdir)
  publishDir "$outputdir", mode:'copy'
  output :
   path("*.vcf.gz"), emit : vcf_gz
   path("*.gz"), emit : all
  script :
    passd = (params.vcf_folder_zip_password!="") ? " -P ${params.vcf_folder_zip_password} " : ""
    """
    unzip $passd $vcfzip
    """
}


process computedstat{
 label 'py3utils'
 memory params.low_memory
  time   params.big_time
  input :
     path(vcf)
     val(vcf_patstatfreq)
     val(vcf_patscoreimp)
  output :
     path("${Ent}")
  script :
    Ent=vcf.baseName+".stat"
    """
    bcftools query -f '%CHROM %REF %ALT %POS %INFO/$vcf_pat %INFO/${vcf_patscoreimp} ${vcf_patstatfreq}\\n' $vcf > $Ent
    """
}

process dostat{
 memory params.low_memory
 input :
    path(all_file)
    val(outputdir)
    val(outputpat)
 publishDir "${outputdir}/", overwrite:true, mode:'copy'
 output :
   tuple path("${fileout}.tex"), emit : tex
   tuple path('*report_freq.pdf'), path("*_report_score.pdf"), emit : support_file
   path("${fileout}*"), emit : all
 script :
  fileout=outputpat+"_report"
  allfile=all_file.join(',')
  """
  stat_vcf_v2.py  --out $fileout --min_score ${params.impute_info_cutoff} --list_files $allfile
  """
}

process clean_vcf {
 label 'py3utils'
 cache 'lenient'
 cpus 5
 input :
   tuple path(vcf), path(fasta),val(maf), val(hwe), val(r2), val(missing), val(outputdir)
 publishDir "${outputdir}", mode:'copy'
 output :
    path(outputfile)
 script :
   cuthwe=hwe > 0 ? " --hwe ${hwe}" : ""
   cutmaf=maf > 0 ? " --maf ${maf}" : ""
   outputfile=vcf.toString().replaceAll(/.vcf.gz/,'_clean.vcf.gz')
   if((cuthwe> 0 || cutmaf>0) & r2>0) {
     """
     vcftools --gzvcf $vcf $cuthwe  $cutmaf --recode --recode-INFO-all --stdout  | ${params.bin_bcftools} view -i '${params.vcf_patscoreimp}>$r2'   | bcftools norm -Ou -m -any  | bcftools norm -Oz -f $fasta  -o $outputfile
     """
   }else if (r2>0){
     """
     ${params.bin_bcftools} view -i '${params.vcf_patscoreimp}>$r2' $vcf  | bcftools norm -Ou -m -any  | bcftools norm -Oz -f $fasta -o $outputfile
     """
   }else {
     """
     vcftools --gzvcf $vcf $cuthwe  $cutmaf --recode --recode-INFO-all --stdout | bcftools norm -Ou -m -any  | bcftools norm -Oz -f $fasta  -o  $outputfile
     """
   }
}

process convert_inplink {
   cpus params.max_cpus
   input :
     tuple path(vcf),val(Ent),val(outputdir), val(nbfile) 
  publishDir "${outputdir}", mode:'copy'
  output :
      tuple path("${Ent}.bed"), path("${Ent}.bim"), path("${Ent}.fam"), emit: bed
      path("${Ent}.bim"), emit: bim
  script :
    if(nbfile>1){
      Ent=vcf.baseName
    }
    """
     plink --vcf $vcf   --vcf-idspace-to _ --double-id --allow-extra-chr --make-bed --out ${Ent} --keep-allele-order 
    """
}

process GetRsDup{
    input :
      path(bim)
    output :
      tuple path(outdel),path(out)
    script :
      lbim=bim.join(",")
      out="snpfile_red.rs"
      outdel="snpfile_del.rs"
      """
      search_dup_bim.py $lbim $out $outdel
      """
}

process TransformRsDup{
    input :
     tuple path(delrange),path(rstochange), path(bed),path(bim),path(fam)
    output :
       tuple path("${newheader}.bed"),path("${newheader}.bim"),path("${newheader}.fam"), emit : plk
       val("$newheader"), emit : plkHead
    script :
       header=bed.baseName
       newheader=header+"_rsqc"
       """
       cp $fam ${header}_save.fam
       cp $bed ${header}_save.bed
       replacers_forbim.py ${bim} $rstochange $header"_save.bim"
       plink -bfile ${header}_save --keep-allele-order --make-bed -out $newheader --exclude range $delrange --allow-extra-chr
       rm ${header}_save*
       """
}


process AddedCM{                                                                
    cpus params.max_cpus
    memory params.plink_mem_req                                                 
    time   params.big_time                                                      
    input :                                                                     
       tuple path(map),path(bedi),path(bimi),path(fami)                         
    output :                                                                    
       tuple path(bedf),path(bimf),path(famf), emit : plk                       
       val("$headerout"),  emit :plkHead                                               
    script :                                                                    
       headerouti=bedi.baseName                                                    
       headerout=headerouti+"_map"                                                    
       bedf=headerout+".bed"                                                       
       bimf=headerout+".bim"                                                       
       famf=headerout+".fam"                                                       
       cm_shap=headerout+".shape"                                                  
       """                                                                      
       chro=`head $bimi|awk '{print \$1}'|uniq`                                 
       sed '1d' $map|awk -v chro=\$chro '{if(chro==\$1)print \$2"\\t"\$3"\\t"\$4}' >> $cm_shap
       awk '{print \$2}' $headerouti".bim" | sort | uniq -d > duplicated_snps.snplist
       plink --bfile $headerouti --exclude duplicated_snps.snplist --make-bed --keep-allele-order --out $headerouti"_tmp" --allow-extra-chr
       plink --bfile $headerouti"_tmp" --list-duplicate-vars ids-only suppress-first --allow-extra-chr
       plink --bfile $headerouti"_tmp" --keep-allele-order --cm-map $cm_shap \$chro   --threads ${params.max_cpus} --make-bed --out $headerout  --exclude plink.dupvar --allow-extra-chr
       """                                                                      
                                                                                
} 
/*
hgrefi=Channel.fromPath(params.reffasta, checkIfExists:true)

process checkfasta{
  cpus params.max_cpus
  label 'py3utils'
  errorStrategy { task.exitStatus == 1 ? 'retry' : 'terminate' }
  maxRetries 1
  publishDir "${params.output_dir}/fastaclean", overwrite:true, mode:'copy'
  input :
     path(fasta) from hgrefi
  output :
     tuple path("$fasta2"), path("${fasta2}.fai") into ref_ch
  script :
    fasta2=fasta.baseName+'_clean.fa.gz'
    """
    if [ "${task.attempt}" -eq "1" ]
    then
    samtools faidx $fasta
    cp $fasta $fasta2
    mv $fasta".fai"  $fasta2".fai"
    else
    zcat $fasta | awk '{print \$1}' | bgzip -@ ${params.max_cpus} -c > $fasta2
    samtools faidx $fasta2
    fi
    """
}

/*
process checkfasta {
  cpus params.max_cpus
  label 'py3utils'
  publishDir "${params.output_dir}/fastaclean", overwrite:true, mode:'copy'
  input :
    path(fasta) from hgrefi
  output :
     tuple path("$fasta2"), path("${fasta2}.fai") into ref_ch
  script :
    fasta2=fasta.baseName+'_clean.fa.gz'
  """
    zcat $fasta | awk '{print \$1}' | bgzip -@ ${params.max_cpus} -c > $fasta2
    samtools faidx $fasta2
  """
}
*/
//
//list_vcf=ref_ch.combine(list_vcf)
//
//if(params.min_scoreinfo>0){
// if(params.lim_used_ram==0){
//  process formatvcfscore{
//   label 'py3utils'
//   cpus params.max_cpus
//   memory { strmem(params.bcftools_mem_req) + 5.GB * (task.attempt -1) }
//   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
//   maxRetries 10
//   time   params.big_time
//   input :
//      tuple path(ref), path(refidx),path(vcf) from list_vcf
//   publishDir "${params.output_dir}/bcf_filter", overwrite:true, mode:'copy', pattern: "*.bcf"
//   output :
//      tuple path("${Ent}.bed"), path("${Ent}.bim"), path("${Ent}.fam") into listchroplink_init
//      path ("${Ent}.bim") into listbimplink_init
//      tuple val(Ent), path("${Ent}.bcf") into bcf_ch
//   script :
//     Ent=vcf.baseName.replaceAll(/.vcf$/, '')+'_filter'
//     threadn=params.max_cpus/ 5
//     threadn = threadn.round()
//     """
//     ln -s $ref tmp.fasta.gz
//     cp $refidx tmp.fasta.gz.fai
//     bcftools view -Ou -i '${params.score_imp}>${params.min_scoreinfo}' $vcf --threads $threadn | bcftools norm -Ou -m -any --threads  $threadn | bcftools norm -Ou -f tmp.fasta.gz --threads $threadn  |bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' --threads  $threadn > $Ent".bcf"
//   cat $Ent".bcf" |plink --bcf /dev/stdin
//     $plink_kao
//     --vcf-idspace-to _
//     --const-fid
//     --allow-extra-chr 0
//     --make-bed
//     --out ${Ent}
//     cp ${Ent}.bim ${Ent}.save.bim
//     awk \'{if(\$2==\".\"){\$2=\$1\":\"\$4\"_\"\$5\"_\"\$6};\$2=substr(\$2,1,20);print \$0}\' ${Ent}.save.bim > ${Ent}.bim
//     """
//   }
//  }else {
//  process formatvcfscore_lowram{
//   label 'py3utils'
//   cpus params.max_cpus
//   memory { strmem(params.bcftools_mem_req) + 5.GB * (task.attempt -1) }
//   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
//   maxRetries 10
//   time   params.big_time
//   input :
//     set file(ref), file(refidx),file(vcf) from list_vcf
//   publishDir "${params.output_dir}/bcf_filter", overwrite:true, mode:'copy', pattern: "*.bcf"
//   output :
//     set file("${Ent}.bed"), file("${Ent}.bim"), file("${Ent}.fam") into listchroplink_init
//     file("${Ent}.bim") into listbimplink_init
//     set val(Ent), file("${Ent}.bcf") into bcf_ch
//   script :
//    Ent=vcf.baseName.replaceAll(/.vcf$/, '')+'_filter'
//    //threadn=params.max_cpus/ 5
//    threadn = params.max_cpus
//    """
//    ln -s $ref tmp.fasta.gz
//    cp $refidx tmp.fasta.gz.fai
//    bcftools view -Ou -i '${params.score_imp}>${params.min_scoreinfo}' $vcf --threads $threadn > vcftmp1
//   cat vcftmp1 | bcftools norm -Ou -m -any --threads  $threadn > vcftmp2
//   rm -f vcftmp1
//   cat vcftmp2 | bcftools norm -Ou -f tmp.fasta.gz --threads $threadn > vcftmp3
//   rm -f vcftmp2
//   cat vcftmp3 |bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' --threads  $threadn > $Ent".bcf"
//   rm vcftmp3
//   cat $Ent".bcf" |plink --bcf /dev/stdin
//    $plink_kao
//    --vcf-idspace-to _
//    --const-fid
//    --allow-extra-chr 0
//    --make-bed
//    --out ${Ent}
//    cp ${Ent}.bim ${Ent}.save.bim
//    awk \'{if(\$2==\".\"){\$2=\$1\":\"\$4\"_\"\$5\"_\"\$6};\$2=substr(\$2,1,20);print \$0}\' ${Ent}.save.bim > ${Ent}.bim
//    """
//  }
// }
//
//process convertbcf_invcf{
//  label 'py3utils'
//  time   params.big_time
//  memory params.bcftools_mem_req
//  publishDir "${params.output_dir}/vcf_filter", overwrite:true, mode:'copy', pattern: "*.vcf.gz"
//  input :
//      set val(Ent), file(bcf) from bcf_ch
//   output :
//      val(vcf)
//   script :
//     vcf=Ent+".vcf.gz"
//     """
//     bcftools convert $bcf -O z > $vcf
//     """
//}
//
//
//}else{
//process formatvcf{
//  label 'py3utils'
//  cpus params.max_cpus
//  memory params.bcftools_mem_req
//  time   params.big_time
//  input :
//     tuple path(ref),path(refidx),path(vcf) from list_vcf
//  output :
//     set file("${Ent}.bed"), file("${Ent}.bim"), file("${Ent}.fam") into listchroplink_init
//     file("${Ent}.bim") into listbimplink_init
//  script :
//    Ent=vcf.baseName+'_filter'
//    """
//    ln -s $ref tmp.fasta.gz
//    cp $refidx tmp.fasta.gz.fai
//    bcftools view -Ou $vcf | bcftools norm -Ou -m -any | bcftools norm -Ou -f tmp.fasta.gz |bcftools annotate -Ob -x ID -I +'%CHROM:%POS:%REF:%ALT' |
//  plink --bcf /dev/stdin
//    $plink_kao
//    --vcf-idspace-to _
//    --const-fid
//    --allow-extra-chr 0
//    --make-bed
//    --out ${Ent}
//
//
//    cp ${Ent}.bim ${Ent}.save.bim
//    awk \'{if(\$2==\".\"){\$2=\$1\":\"\$4\"_\"\$5\"_\"\$6};\$2=substr(\$2,1,20);print \$0}\' ${Ent}.save.bim > ${Ent}.bim
//    """
//}
//}
//
//if(params.file_ref_gzip!=""){
// rs_infogz=Channel.fromPath(params.file_ref_gzip,checkIfExists:true)
//listchroplink_init_ch2=listchroplink_init.combine(rs_infogz)
process extractrsname{
  memory params.high_memory
  cpus params.cpus_
  input :
    tuple path(bed), path(bim), path(fam), path(rsinfo), val(outputdir)
  publishDir "$outputdir",  mode:'copy', pattern: "*_pos"
  output :
    tuple path("${out}.bed"), path("${out}.bim"), path("${out}.fam"), emit: plk
    path("${out}.bim"), emit:  bim
    path("$outrs*"), emit : log
  script :
    plk=bed.baseName
    out=plk+'_rsupdate'
    outrs=bed.baseName+'_pos'
    extract=params.deleted_notref=='F' ? "" : " --extract range keep"
    """
   zcat $rsinfo | extractrsid_bypos.py --bim $bim --out_file $outrs --ref_file stdin --chro_ps ${params.poshead_chro_inforef} --bp_ps ${params.poshead_bp_inforef} --rs_ps ${params.poshead_rs_inforef} --a1_ps ${params.poshead_a1_inforef}  --a2_ps ${params.poshead_a2_inforef}
   awk '{print \$1"\t"\$2"\t"\$2"\t"\$5}' $outrs > keep
   plink --keep-allele-order --noweb $extract --bfile $plk --make-bed --out $out --update-name $outrs".rs" -maf 0.0000000000000000001 --threads ${params.max_cpus}
    """
}
//
//
//}else{
//listbimplink=listbimplink_init
//listchroplink=listchroplink_init
//}
//
//
//
//if(params.genetic_maps!=""){
//GMMap=Channel.fromPath(params.genetic_maps, checkIfExists:true)
//listchroplinkrsmap=GMMap.combine(listchroplinkrs)
//process AddedCM{
//    cpus params.max_cpus
//    memory params.plink_mem_req
//    time   params.big_time
//    input :
//       set file(map),file(bedi),file(bimi),file(fami) from listchroplinkrsmap
//    output :
//       set file(bedf),file(bimf),file(famf) into listchroplinkrsf
//       val(header) into plinkheadf
//    script :
//       headeri=bedi.baseName
//       header=headeri+"_map"
//       bedf=header+".bed"
//       bimf=header+".bim"
//       famf=header+".fam"
//       cm_shap=header+".shape"
//       """
//       chro=`head $bimi|awk '{print \$1}'|uniq`
//       sed '1d' $map|awk -v chro=\$chro '{if(chro==\$1)print \$2"\\t"\$3"\\t"\$4}' >> $cm_shap
//       awk '{print \$2}' $headeri".bim" | sort | uniq -d > duplicated_snps.snplist
//       plink --bfile $headeri --exclude duplicated_snps.snplist --make-bed $plink_kao --out $headeri"_tmp"
//       plink --bfile $headeri"_tmp" --list-duplicate-vars ids-only suppress-first
//       plink --bfile $headeri"_tmp" $plink_kao --cm-map $cm_shap \$chro   --threads ${params.max_cpus} --make-bed --out $header  --exclude plink.dupvar
//       rm $headeri"_tmp"*
//       """
//
//}
//
//}else{
//plinkheadf=plinkhead
//listchroplinkrsf=listchroplinkrs
//}
//listplinkinf=listchroplinkrsf.collect()
//headplinkf=plinkheadf.collect()
//
process MergePlink{
  cpus params.max_cpus
  memory params.high_memory
  time   params.big_time
  input :
       path(lplk)
       val(outputdir)
       val(outputpat)
  publishDir "${outputdir}/", mode:'copy'
  output :
    tuple path("${outputpat}.bed"), file("${outputpat}.bim"),file("${outputpat}.fam") 
  script :
       hplk2=lplk.join(',')
       """
       echo $hplk2 | awk -F',' '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' | sed 's/\\.[^.]*\$//'  | sort |uniq |sed '1d'> fileplk
       hplkFirst=`echo $hplk2 | awk -F',' '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' | sed 's/\\.[^.]*\$//'  | sort |uniq |head -1`
       plink --bfile \$hplkFirst --keep-allele-order --threads ${params.max_cpus} --merge-list fileplk --make-bed --out $outputpat --allow-extra-chr
       """
}
//
// process clean_plink {
//  cpus params.max_cpus
//  memory params.plink_mem_req
//  time   params.big_time
//  input:
//       tuple path(bed), path(bim), path(fam) from plinkformatf_merge
//  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
//  output :
//     tuple path("${output}.bed"), path("${output}.bim"),path("${output}.fam") into plinkformatf
//  script :
//  plk=bed.baseName
//  output = params.output_pat+'_clean'
//  sexinfo = "--allow-no-sex"
//  """
//  if [ $params.cut_mind -eq 0 ]
//  then
//  bfile=$plk
//  else
//  plink $plink_kao --bfile $plk $sexinfo --mind $params.cut_mind --make-bed --out temp2
//  bfile=temp2
//  fi
//  if [ $params.cut_geno -eq 0 ]
//  then
//  bfile2=\$bfile
//  else
//  plink $plink_kao --bfile \$bfile   $sexinfo --geno $params.cut_geno --make-bed --out temp3
//  bfile2=temp3
//  /bin/rm -f temp2.{bim,fam,bed}
//  fi
//  if [ $params.cut_maf -eq 0 ]
//  then
//  bfile3=\$bfile2
//  else
//  plink $plink_kao --bfile \$bfile2 $sexinfo --maf $params.cut_maf --make-bed --out temp4
//  /bin/rm -f temp3.{bim,fam,bed}
//  /bin/rm -f temp2.{bim,fam,bed}
//  bfile3=temp4
//  fi
//  if [ $params.cut_hwe -eq 0 ]
//  then
//  bfile4=\$bfile3
//  else
//  plink $plink_kao --bfile \$bfile3 $sexinfo --hwe $params.cut_hwe --make-bed  --out temp5
//  /bin/rm -f temp4.{bim,fam,bed}
//  /bin/rm -f temp3.{bim,fam,bed}
//  /bin/rm -f temp2.{bim,fam,bed}
//  bfile4=temp5
//  fi
//  cp \$bfile4".bed" $output".bed"
//  cp \$bfile4".bim" $output".bim"
//  cp \$bfile4".fam" $output".fam"
//  rm -f temp*.{bim,fam,bed}
//  """
//
// }
//
//if(params.data!=''){
//
//data_ch=channel.fromPath(params.data, checkIfExists:true)
process plink_update_name{
  label 'R'
  cpus params.max_cpus
  memory params.high_memory
 input :
   tuple path(bed), path(bim), path(fam)
   path(data)
   val(outputdir)
   val(outputpat)
  publishDir "${outputdir}/", mode:'copy'
 output :
     tuple path("${outputpat}.bed"), path("${outputpat}.bim"),path("${outputpat}.fam")
     path("${out}.log")
 script :
     out=params.output_pat+"_idupdate"
     bfile=bed.baseName
     """
     change_updateidplink.r $fam $data update_id
     plink -bfile $bfile -make-bed $plink_kao --update-ids update_id -out $out --threads ${params.max_cpus} --allow-extra-chr
     """
}
//
//}
