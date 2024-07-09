// clean rs
process extract_bp_fromvcf{
 input :
      path(filebed)
      tuple path(file_annot), path(file_annotindex)
      val(outputdir)
      val(outputname)
 publishDir "${outputdir}/",  mode:'copy'
 output :
     path("${outputname}.vcf")
 script :
   """
   ${params.bin_bcftools} view $file_annot  -R $filebed > $outputname".vcf"
   """
}

process extract_rs_fromvcf{
 input :
      path(filebed)
      tuple path(file_annot), path(file_annotindex)
      val(outputdir)
      val(outputname)
 publishDir "${outputdir}/",  mode:'copy'
 output :
     path("${outputname}.vcf")
 script :
   """
   ${params.bin_bcftools} view  -i'ID=@snplist.txt' $file_annot > ${outputname}".vcf"
   """
}

process bcftools_index_vcf{
 label 'py3utils'
 cpus params.max_plink_cores
 input :
   path(filegz)
 output :
   tuple path(filegz), path("${filegz}.csi")
 output :
   """
   ${params.bin_bcftools} index  $filegz --threads ${params.max_plink_cores}
   """
}

process checkfixref{
  label 'py3utils'
  input :
    path(vcf)
    tuple path(hg), path(index)
  publishDir "${params.output_dir}/check/Bcftools", overwrite:true, mode:'copy'
  output :
    path("${params.output}.checkbcf*")
  script :
    """
    ${params.bin_bcftools} +fixref $vcf -- -f $hg 1> ${params.output}".checkbcf.out" 2> ${params.output}".checkbcf.err"
    """
}
 process list_chro {                                                            
  input :                                                                       
      path(vcf)
  output :                                                                      
     stdout                                                                     
  script :                                                                      
      """                                                                       
      zcat $vcf | grep -v "#" |awk '{print \$1}' |uniq|sort|uniq                                     
      """                                                                       
 }   
process splitvcf {
  input :
    path(vcf)
    val(outputdir)
    val(output)
  publishDir "${outputdir}/", mode:'copy'
  output :
   path("output*.vcf.gz")
  script :
    """
    zcat $vcf | split_chrovcf.py $output
    """
}
