
process clean_phenofile {
  label 'R'
  input :
    path(data)
    tuple path(bed), path(bim), path(fam)
  output :
    path("${out}.dup"), emit : dup
    path("${out}.corname"), emit : correspond
    tuple path("${out}_phenoind_notfound.csv"), path("${out}_fam_notfound.csv"), emit : errorind
    tuple path("${out}.update"), emit : pheno_i
  script :
   out="clean_dup"
   bfile=bed.baseName
   """
   extract_duplicate.r --data $data --bfile $bfile --out $out --col_fidid ${params.col_fidid}  --col_newfidid ${params.newcol_fidid}
   """ 
}

process extract_ind_plink{
  input :
    tuple path(bed), path(bim), path(fam)
    path(keep)
  output :
    tuple path("${out}.bed"), path("${out}.bim"), path("${out}.fam")
  script :
   bfile=bed.baseName
   out=bfile+'_clean'
  """
  plink -bfile $bfile --keep  $keep --make-bed -out $out --keep-allele-order --maf 0.01   --geno 0.01 --autosome
  """
}

process plink_indep{
  input :
    tuple path(bed), path(bim), path(fam), val(outputdir), val(out)
  output :
   tuple path("${out}.prune.in"), emit : list_pos
   tuple path("${out}.bed"), path("${out}.bim"), path("${out}.fam"), emit : plk
 script :
   bfile=bed.baseName                                                           
  """
  plink -bfile $bfile  -out $out --indep-pairwise 100 20 0.2 
  plink -bfile $bfile --extract ${out}.prune.in --keep-allele-order -make-bed -out $out
  """
}

process computed_relatdness {
   input :
     tuple path(bed), path(bim), path(fam)
     tuple val(isfileind),path(fileind)
     val(full)
     val(minrel)
     val(out)
   output :
    path("${out}.genome")
  script :
   bfile=bed.baseName
   full=""
   if(full.toString()=='T')full=" full" 
   cmdind=(isfileind==1) ?  "--keep $fileind" : ""
   minrel=(minrel>0) ? " --min 0.7 " : ""
   """
   plink -bfile $bfile -out $out $cmdind $minrel --genome $full
   """
}

process compute_missing{
  input :
     tuple path(bed), path(bim), path(fam)
     tuple val(isfileind),path(fileind)
     val(out)
  output :
    path("${out}.imiss"),emit : ind
    path("${out}.lmiss"),emit : snp
  script :
    bfile=bed.baseName
    cmdind=(isfileind==1) ?  "--keep $fileind" : ""
    """
    plink -bfile $bfile --out $out $cmdind --missing
    """
}

process check_rel {
  label 'R'
  input :
   path(corrname)
   path(datai)
   path(relall) 
   path(reldup) 
   path(miss)
   val(minpihat)
   val(out)
   val(outputdir)
  publishDir "${outputdir}/", overwrite:true, mode:'copy'
  output :
    path("${out}.pheno"), emit : pheno
    path("${out}.indtokeep"), emit : indtokeep
    path("${out}.id_update"), emit : update_id
    path("$out*"), emit : all
  script :
  """
 check_rel.r --corname $corrname --rel_all $relall --rel_dup $reldup  --missing_dup $miss --min_pihat $minpihat --out $out  --data_i $datai
 """

}

process plink_updatename {
 input :
   tuple path(bed), path(bim), path(fam)                                      
   path(fileupdate)
   val(out)
   val(outputdir)
 publishDir "${outputdir}/", overwrite:true, mode:'copy'
 output :
   tuple path("${out}.bed"), path("${out}.bim"), path("${out}.fam")
 script :
   bfile = bed.baseName
   """
   awk '{print \$3"\t"\$4}' $fileupdate  > keep
   plink -bfile $bfile --update-ids $fileupdate --make-bed -out $out --keep-allele-order --keep keep
   """

}
