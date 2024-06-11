
process updateplk_rsname{
  label 'py3utils'
  memory params.plink_mem_req
  cpus params.max_plink_cores
  input :
    tuple path(rsinfo), path(rsinfo_csi)
    tuple path(bed), path(bim), path(fam)
    val(outputdir)
  publishDir "${outputdir}", overwrite:true, mode:'copy'
  output :
    path(outrs), emit : rs_info
    tuple path("${out}.bed"),path("${out}.bim"),path("${out}.fam"), emit : plink
  script :
    plk=bed.baseName
    outrs=plk+"_updaters"
    out=plk+"_updaters"
    extract=params.rsinfovcf_deleted_notref=='F' ? "" : " --extract range keep"
    """
    zcat $rsinfo| extractrsid_bypos.py --bim $bim --out_file $outrs --ref_file stdin --chro_ps ${params.rsinfovcf_poshead_chro_inforef} --bp_ps ${params.rsinfovcf_poshead_bp_inforef} --rs_ps ${params.rsinfovcf_poshead_rs_inforef} --a1_ps ${params.rsinfovcf_poshead_a1_inforef}  --a2_ps ${params.rsinfovcf_poshead_a2_inforef}
    awk '{if(\$1=="X" || \$1=="chrX"){\$1=23};print \$1"\t"\$2"\t"\$2"\t"\$5}' $outrs > keep
    plink --keep-allele-order $extract --bfile $plk --make-bed --out $out --update-name $outrs".rs" -maf 0.0000000000000000001 --threads ${params.max_plink_cores}
    """
}

process updateplk_rsname_norsinfo{
  label 'py3utils'
  memory params.plink_mem_req
  cpus params.max_plink_cores
  input :
    tuple path(bed), path(bim), path(fam)
    val(outputdir)
  publishDir "${outputdir}", overwrite:true, mode:'copy'
  output :
    path("${out}.bed"),path("${out}.bim"),path("${out}.fam"), emit : plink
  script :
    plk=bed.baseName
    outrs=plk+"_updaters"
    out=plk+"_updaters"
    extract=params.rsinfovcf_deleted_notref=='F' ? "" : " --extract range keep"
    """
    zcat $rsinfo| extractrsid_bypos.py --bim $bim --out_file $outrs --ref_file stdin --chro_ps ${params.poshead_chro_inforef} --bp_ps ${params.poshead_bp_inforef} --rs_ps ${params.poshead_rs_inforef} --a1_ps ${params.poshead_a1_inforef}  --a2_ps ${params.poshead_a2_inforef}
    awk '{if(\$1=="X" || \$1=="chrX"){\$1=23};print \$1"\t"\$2"\t"\$2"\t"\$5}' $outrs > keep
    plink --keep-allele-order $extract --bfile $plk --make-bed --out $out --update-name $outrs".rs" -maf 0.0000000000000000001 --threads ${params.max_plink_cores}
    """
}


process deletedmultianddel{
   label 'R'
   memory params.plink_mem_req
   cpus params.max_plink_cores
   input :
    tuple path(bed), path(bim), path(fam)
   output :
    tuple path("${out}.bed"),path("${out}.bim"),path("${out}.fam")
   script :
    plk=bed.baseName
    out=plk+"_nomulti"
    """
    biv_selgoodallele.r $bim rstodel ${params.convertinvcf_justagtc}
    plink --keep-allele-order --make-bed --bfile $plk --out $out -maf 0.0000000000000000001 --exclude rstodel --threads ${params.max_plink_cores}
    """
}

process refallele{
   memory params.plink_mem_req
   cpus params.max_plink_cores
   input :
    tuple path(bed), path(bim), path(fam)
    path(infors)
   output :
    tuple path("${out}.bed"),path("${out}.bim"),path("${out}.fam")
  script :
    plk=bed.baseName
    out=plk+"_refal"
    """
    awk '{print \$5"\t"\$6}' $infors > alleref
    plink --bfile $plk --make-bed --out $out --threads ${params.max_plink_cores} --a2-allele alleref
    """
}



process mergevcf{
  label 'py3utils'
  cpus params.max_plink_cores
  input :
   path(allfile)
   val(outputdir)
  publishDir "${outputdir}/", overwrite:true, mode:'copy'
  output :
     path("${out}.vcf.gz")
  script :
    fnames = allfile.join(" ")
    out="${params.output}"
    """
    ${params.bin_bcftools} concat -Oz -o ${out}.vcf.gz --threads ${params.max_plink_cores} $fnames
    """
}

process checkVCF{
  label 'py3utils'
  input :
    path(vcf)
    tuple path(hg), path(index)
    val(outputdir)
  publishDir "${outputdir}", overwrite:true, mode:'copy'
  output :
    path("${out}*")
  script :
    out="${params.output}_check"
    """
    checkVCF.py -r $hg -o $out $vcf
    """
}


/*
process extractrsname{
  label 'py3utils'
  memory params.plink_mem_req
  cpus params.max_plink_cores
  input :
    tuple path(rsinfo), path(rsinfo_csi)
    tuple path(bed), path(bim), path(fam)
  publishDir "${params.output_dir}/rsclean", overwrite:true, mode:'copy'
  output :
    tuple path("${out}.bed"),path("${out}.bim"),path("${out}.fam"), emit: plk
    path(outrs), emit : rs_info
  script :
    plk=bed.baseName
    outrs=plk+"_updaters"
    out=plk+"_updaters"
    extract=params.rsinfovcf_deleted_notref=='F' ? "" : " --extract range keep"
    """
    zcat $rsinfo| extractrsid_bypos.py --bim $bim --out_file $outrs --ref_file stdin --chro_ps ${params.rsinfovcf_poshead_chro_inforef} --bp_ps ${params.rsinfovcf_poshead_bp_inforef} --rs_ps ${params.rsinfovcf_poshead_rs_inforef} --a1_ps ${params.rsinfovcf_poshead_a1_inforef}  --a2_ps ${params.rsinfovcf_poshead_a2_inforef}
    awk '{if(\$1=="X" || \$1=="chrX"){\$1=23};print \$1"\t"\$2"\t"\$2"\t"\$5}' $outrs > keep
    plink --keep-allele-order $extract --bfile $plk --make-bed --out $out --update-name $outrs".rs" -maf 0.0000000000000000001 --threads ${params.max_plink_cores}
    """
}
*/

process convertInVcf {
   label 'py3utils'
   cache 'lenient'
   memory params.low_memory
   cpus params.max_plink_cores
   time params.big_time
   input :
    tuple path(bed), path(bim), path(fam), path(gz_info), path(gz_info_csi), path(fast), path(fastaindex),val(chro), val(outputdir)
   publishDir "$outputdir/", overwrite:true, mode:'copy'
   output :
    path("${out}.vcf.gz"), emit : vcf
    path("${out}.rep"), emit : rep
   script:
     base=bed.baseName
     parmchro=""
     if(chro!=""){
       parmchro="--chr $chro"
       out="${params.output}_${chro}"
     }else{
       out="${params.output}"

    }
     """
     mkdir -p ${params.tmpdir}
     plink2  --bfile ${base}  --recode vcf-iid bgz --out $out --keep-allele-order --snps-only --threads ${params.max_plink_cores}  $parmchro
     ${params.bin_bcftools} view ${out}.vcf.gz | bcftools sort - -O z -T ${params.tmpdir} > ${out}_tmp.vcf.gz
     rm -f ${out}.vcf.gz
     ${params.bin_bcftools} +fixref ${out}_tmp.vcf.gz -Oz -o ${out}.vcf.gz -- -f $fast -m flip -d &> $out".rep"
     """
 }

process checkfixref{
  label 'py3utils'
  input :
    path(vcf)
    tuple path(hg), path(index)
    val("outputdir")
  publishDir "${params.output_dir}/check/Bcftools", overwrite:true, mode:'copy'
  output :
    path("${params.output}.checkbcf*")
  script :
    """
    ${params.bin_bcftools} +fixref $vcf -- -f $hg 1> ${params.output}".checkbcf.out" 2> ${params.output}".checkbcf.err"
    """
}
