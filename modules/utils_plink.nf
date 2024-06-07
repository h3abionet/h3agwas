
 process getfrequency{
 input :
     tuple path(bed), path(bim), path(fam)
 publishDir "${params.output_dir}/checkmichigan/frq", overwrite:true, mode:'copy'
 output :
     tuple path(bim),path("${bfile}.frq")
 script :
    bfile=bed.baseName
    """
    plink --freq -bfile $bfile  -out $bfile
    """
 }

 process list_chro {
  input :
      tuple path(bed), path(bim), path(fam)                                      
  publishDir "${outputdir}/", overwrite:true, mode:'copy'
  output :
     stdout
  script :
      """
      awk '{print \$1}' $bim|uniq|sort|uniq                                                
      """
 }
