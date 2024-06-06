
 process clean_michigan{
  input :
    tuple path(exc), path(chro), path(pos), path(strand), path(Force)
    tuple path(bed), path(bim), path(fam)
  output :
    tuple path("${out}.bed"), path("${out}.bim"), path("${out}.fam")
  script :
  base=bed.baseName
  out=base+"_updated"
  """
  plink --bfile $base --exclude $exc --make-bed --out TEMP1
  plink --bfile TEMP1 --update-map $chro --update-chr --make-bed --out TEMP2
  plink --bfile TEMP2 --update-map $chro --make-bed --out TEMP3
  plink --bfile TEMP3 --flip $strand --make-bed --out TEMP4
  plink --bfile TEMP4 --reference-allele $Force --make-bed --out $out
  rm TEMP*
  """
 }


 process michigan_qc {
  input :
      tuple path(bim), path(frq)
      path(binmich)
      path(datamich)
  publishDir "${params.output_dir}/checkmichigan/output", overwrite:true, mode:'copy'
  output :
      tuple path("Exclude-$base-HRC.txt"), path("Chromosome-*-HRC.txt"), path("Position-$base-HRC.txt"), path("Strand*HRC.txt"), path("Force-*-HRC.txt"), emit : res
      path("LOG-$base*"), emit : log
  script :
    base=bim.baseName
    """
    $binmich -b $bim -f $frq -r $datamich -h
    """
 }

 process dl_dataref_michigan{
    publishDir "${params.output_dir}/checkmichigan/data", overwrite:true, mode:'copy'
    output :
      path("*.tab")
    script :
    """
    wget ${params.ftp_dataref_michigan}
    gunzip `basename ${params.ftp_dataref_michigan}`
    """
 }
