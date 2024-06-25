#!/usr/bin/env nextflow

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
  output :
     path("${Ent}")
  script :
    Ent=vcf.baseName+".stat"
    """
    bcftools query -f '%CHROM %REF %ALT %POS %INFO/$params.vcf_pat %INFO/${params.vcf_patscoreimp} ${params.vcf_patstatfreq}\\n' $vcf > $Ent
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
  stat_vcf_v2.py  --out $fileout --min_score ${params.vcf_minscoreimp} --list_files $allfile
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
   cpus params.max_cpu 
   input : 
     path(vcf)
     val(Ent)
     val("outputdir")
  publishDir "${outputdir}", mode:'copy'
  output :
      tuple path("${Ent}.bed"), path("${Ent}.bim"), path("${Ent}.fam") 
  script :
    """
     plink --vcf $vcf   --vcf-idspace-to _ --double-id --allow-extra-chr --make-bed --out ${Ent} --keep-allele-order
    """
}
/*
hgrefi=Channel.fromPath(params.reffasta, checkIfExists:true)

process checkfasta{
  cpus params.max_plink_cores
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
    zcat $fasta | awk '{print \$1}' | bgzip -@ ${params.max_plink_cores} -c > $fasta2
    samtools faidx $fasta2
    fi
    """
}


/*
process checkfasta {
  cpus params.max_plink_cores
  label 'py3utils'
  publishDir "${params.output_dir}/fastaclean", overwrite:true, mode:'copy'
  input :
    path(fasta) from hgrefi
  output :
     tuple path("$fasta2"), path("${fasta2}.fai") into ref_ch
  script :
    fasta2=fasta.baseName+'_clean.fa.gz'
  """
    zcat $fasta | awk '{print \$1}' | bgzip -@ ${params.max_plink_cores} -c > $fasta2
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
//   cpus params.max_plink_cores
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
//     threadn=params.max_plink_cores/ 5
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
//   cpus params.max_plink_cores
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
//    //threadn=params.max_plink_cores/ 5
//    threadn = params.max_plink_cores
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
//  cpus params.max_plink_cores
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
// process extractrsname{
//  memory params.plink_mem_req
//  cpus params.max_plink_cores
//  input :
//    set file(bed), file(bim), file(fam), file(rsinfo) from listchroplink_init_ch2
//  publishDir "${params.output_dir}/rs/extractrsname", overwrite:true, mode:'copy', pattern: "*_pos"
//  output :
//    set path("${out}.bed"), path("${out}.bim"), path("${out}.fam") into listchroplink
//    path("${out}.bim") into listbimplink
//    path("$outrs*")
//  script :
//    plk=bed.baseName
//    out=plk+'_rsupdate'
//    outrs=bed.baseName+'_pos'
//    extract=params.deleted_notref=='F' ? "" : " --extract range keep"
//    """
//   zcat $rsinfo | extractrsid_bypos.py --bim $bim --out_file $outrs --ref_file stdin --chro_ps ${params.poshead_chro_inforef} --bp_ps ${params.poshead_bp_inforef} --rs_ps ${params.poshead_rs_inforef} --a1_ps ${params.poshead_a1_inforef}  --a2_ps ${params.poshead_a2_inforef}
//   awk '{print \$1"\t"\$2"\t"\$2"\t"\$5}' $outrs > keep
//   plink --keep-allele-order --noweb $extract --bfile $plk --make-bed --out $out --update-name $outrs".rs" -maf 0.0000000000000000001 --threads ${params.max_plink_cores}
//    """
// }
//
//
//}else{
//listbimplink=listbimplink_init
//listchroplink=listchroplink_init
//}
//
//
//
//bimmerg=listbimplink.collect()
//process GetRsDup{
//    memory params.other_mem_req
//    input :
//      file(bim) from bimmerg
//    publishDir "${params.output_dir}/rs/getrsdup", overwrite:true, mode:'copy'
//    output :
//      set file(outdel),file(out) into duplicat
//    script :
//      lbim=bim.join(",")
//      out="snpfile_red.rs"
//      outdel="snpfile_del.rs"
//      """
//      search_dup_bim.py $lbim $out $outdel
//      """
//}
//listchroplinkinf=duplicat.combine(listchroplink)
//process TransformRsDup{
//    input :
//     set file(delrange),file(rstochange), file(bed),file(bim),file(fam) from listchroplinkinf
//    publishDir "${params.output_dir}/rs/transformrsdup", overwrite:true, mode:'copy'
//    output :
//       set file("${newheader}.bed"),file("${newheader}.bim"),file("${newheader}.fam") into listchroplinkrs
//       val(newheader) into plinkhead
//    script :
//       header=bed.baseName
//       newheader=header+"_rsqc"
//       """
//       cp $bed ${header}.save.bed
//       cp $bim ${header}.save.bim
//       cp $fam ${header}.save.fam
//       cp ${header}.save.bim  ${bim}.save
//       replacers_forbim.py ${bim}.save $rstochange ${header}.save.bim
//       plink -bfile ${header}.save $plink_kao --make-bed -out $newheader --exclude range $delrange
//       rm  ${header}.save.*
//       """
//}
//
//
//if(params.genetic_maps!=""){
//GMMap=Channel.fromPath(params.genetic_maps, checkIfExists:true)
//listchroplinkrsmap=GMMap.combine(listchroplinkrs)
//process AddedCM{
//    cpus params.max_plink_cores
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
//       plink --bfile $headeri"_tmp" $plink_kao --cm-map $cm_shap \$chro   --threads ${params.max_plink_cores} --make-bed --out $header  --exclude plink.dupvar
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
//process MergePlink{
//  cpus max_plink_cores_merge
//  memory plink_mem_req_merge
//  time   params.big_time
//  input :
//       path(lplk) from listplinkinf
//       val(hplk) from headplinkf
//  publishDir "${params.output_dir}/nofilter", overwrite:true, mode:'copy'
//  output :
//     set path("${output}.bed"), file("${output}.bim"),file("${output}.fam") into plinkformatf_merge
//  script :
//       output=params.output_pat+"_beforefilter"
//       hplk2=lplk.join(',')
//       """
//       echo $hplk2 | awk -F',' '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' | sed 's/\\.[^.]*\$//'  | sort |uniq |sed '1d'> fileplk
//       hplkFirst=`echo $hplk2 | awk -F',' '{for(Cmt=1;Cmt<=NF;Cmt++)print \$Cmt}' | sed 's/\\.[^.]*\$//'  | sort |uniq |head -1`
//       plink --bfile \$hplkFirst --keep-allele-order --threads ${max_plink_cores_merge} --merge-list fileplk --make-bed --out $output
//       """
//}
//
// process clean_plink {
//  cpus params.max_plink_cores
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
//process update_name{
//  label 'R'
//  cpus params.max_plink_cores
//  memory params.plink_mem_req
// input :
//   tuple path(bed), path(bim), path(fam) from plinkformatf
//   path(data) from data_ch
//  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
// output :
//     tuple path("${params.output_pat}_idupdate.bed"), path("${params.output_pat}_idupdate.bim"),path("${params.output_pat}_idupdate.fam") into plinkformatf_chgid
//     path("${out}.log")
// script :
//     out=params.output_pat+"_idupdate"
//     bfile=bed.baseName
//     """
//     change_updateidplink.r $fam $data update_id
//     plink -bfile $bfile -make-bed $plink_kao --update-ids update_id -out $out
//     """
//}
//
//}
