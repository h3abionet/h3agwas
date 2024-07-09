process checkfasta{
  cpus params.max_cpus
  label 'py3utils'
  errorStrategy { task.exitStatus == 1 ? 'retry' : 'terminate' }
  maxRetries 1
  input :
     path(fasta)
     val(outputdir)
 publishDir "${outputdir}/", overwrite:true, mode:'copy'
  output :
     tuple path("$fasta2"), path("${fasta2}.fai"), emit: fasta_index
     path("$fasta2"), emit: fasta
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


process checkfasta_2{
  cpus params.max_cpus
  label 'py3utils'
  input :
     path(fasta)
     val(outputdir)
 publishDir "${outputdir}/", overwrite:true, mode:'copy'
  output :
     tuple path("$fasta2"), path("${fasta2}.fai"), emit: fasta_index
     path("$fasta2"), emit: fasta
  script :
    fasta2=fasta.baseName+"_clean.fa.gz"
    """
    zcat $fasta | bgzip -@ ${params.max_plink_cores} -c > $fasta2
    samtools faidx $fasta2
    """
}

process MD5_plk {
  input:
     path(plink)
  output:
     path(out)
  script:
       bed = plink[0]
       bim = plink[1]
       fam = plink[2]
       out  = "${plink[0].baseName}.md5"
       template "md5_plk.py"
}


process latex_compilation{                                                                
  label 'latex'                                                                 
  input :                                                                       
    path(tex)
    path(otherfile) 
    val(outputdir)
  publishDir "${outputdir}/", overwrite:true, mode:'copy'               
  output :                                                                      
   path("${out}.pdf")                                                           
  script :                                                                      
     out=tex.baseName                                                           
    """                                                                        
    pdflatex $out >& /dev/null                                                  
    pdflatex $out                                                               
    """                                                                         
}
