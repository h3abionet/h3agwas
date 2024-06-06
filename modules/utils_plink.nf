
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

