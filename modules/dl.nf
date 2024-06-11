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
