process dl_fasta{
 label 'py3utils'
 input :
   val(ref)
   val(link)
   val(out)
   val(outputdir)
 publishDir "${outputdir}/", overwrite:true, mode:'copy'
 output :
   tuple path("$out"), emit : fasta
   tuple path("$out"), path("${out}.fai"), emit : fasta_index
 script :
  """
  wget -O - $link |awk '{print \$1}'|bgzip -c > $out
  samtools index $out
  """
}

process dl_gwascat{
   input :
      val(ftp)
      val(out)
      val(outputdir)
   publishDir "${outputdir}",  mode:'copy'
   output :
      path("$out")
   script :
   """
   wget -c $ftp --no-check-certificate -o $out
   """
}

process formatgwascat_pheno{
   label 'R'
   publishDir "${params.output_dir}/gwascat",  mode:'copy'
   input :
      path(gwascat)
      val()
   output :
       file("${out}*")
   script :
     out=params.output
   """
   format_gwascat_pheno.r --file $gwascat --out $out --chro_head $gwascathead_chr --bp_head $gwascathead_bp --pheno_head ${params.head_pheno_gwascat} --beta_head ${params.head_beta_gwascat} --ci_head ${params.head_ci_gwascat} --p_head ${params.head_pval_gwascat}  --n_head ${params.head_n_gwascat} --freq_head ${params.head_af_gwascat} --rs_head ${params.head_rs_gwascat} --riskall_head ${params.head_riskall_gwascat} --format $format --typeformat $typeformat
   """
}
