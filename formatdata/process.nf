
process updateplk_rsname{
  label 'py3utils'
  memory params.plink_mem_req
  cpus params.max_plink_cores
  input :                                                                       
    tuple path(rsinfo), path(rsinfo_csi)
    tuple path(bed), path(bim), path(fam) 
  publishDir "${params.output_dir}/rsclean", overwrite:true, mode:'copy'
  output :                                                                      
    path(outrs), emit : rs_name
    path("${out}.bed"),path("${out}.bim"),path("${out}.fam") 
  script :
    plk=bed.baseName
    outrs=plk+"_updaters"
    out=plk+"_updaters"
    extract=params.deleted_notref=='F' ? "" : " --extract range keep" 
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
    biv_selgoodallele.r $bim rstodel
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


process checkfasta{                                                             
  cpus params.max_plink_cores                                                   
  label 'py3utils'                                                              
  errorStrategy { task.exitStatus == 1 ? 'retry' : 'terminate' }                
  maxRetries 1                                                                  
  input :                                                                       
     path(fasta) 
 publishDir "${params.output_dir}/fasta/", overwrite:true, mode:'copy'          
  output :                                                                      
     tuple path("$fasta2"), path("${fasta2}.fai") 
  script :                                                                      
    fasta2=fasta.baseName+"_clean.fa.gz"                                        
    """                                                                         
    if [ "${task.attempt}" -eq "1" ]                                            
    then                                                                        
    samtools faidx $fasta                                                       
    cp $fasta $fasta2                                                           
    mv $fasta".fai"  $fasta2".fai"                                              
    else                                                                        
    zcat $fasta | bgzip -@ ${params.max_plink_cores} -c > $fasta2               
    samtools faidx $fasta2
    fi
    """                                                                         
}



process mergevcf{                                                               
  label 'py3utils'                                                              
  cpus params.max_plink_cores                                                   
  input :                                                                       
   path(allfile) 
  publishDir "${params.output_dir}/vcf/", overwrite:true, mode:'copy'           
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
  publishDir "${params.output_dir}/check/CheckVCF", overwrite:true, mode:'copy' 
  output :                                                                      
    path("${out}*")                                                             
  script :                                                                      
    out="${params.output}_check"                                                
    """                                                                         
    checkVCF.py -r $hg -o $out $vcf                                             
    """                                                                         
}  
