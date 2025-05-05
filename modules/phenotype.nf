process extractPheno {
    input:
     path(data)
     val(pheno)
     val(covariates)
     val(gxe)
    output:
     path(phenof)
    script:
     phenof = data.baseName+".clean_pheno"
     all_phenos = covariates.length()>0 ? pheno+","+covariates : pheno
     all_phenos = gxe=='' ? all_phenos :  all_phenos+","+gxe
     """
     extract_pheno.py $data ${all_phenos} $phenof
     """
}

process showPhenoDistrib {
    // not sure difference between container and label
    input:
    path(data)
    val(pheno)
    output:
      path("B050*")
    script:
      """
      phe_distrib.py --pheno ${params.pheno} $data B050
      """
}

process add_pcs{
  label 'R'
  cpus params.max_cpus
  input :
      tuple path(eigenval), path(eigenvec), path(phenofile), val(npcs), val(covariates), val(covariates_bin)
  publishDir "${params.output_dir}/pheno/pca/", mode:'copy'
  output :
    path(newfilepheno), emit :data
    env pcs, emit : pcs_covar
    env cov_type, emit :covar_type
  script :
     head=phenofile.baseName
     newfilepheno=head+"_"+npcs.toString()+".newpheno"
     covar2=(covariates=="")? "" : " --covar $covariates "

     """
     addpheno_pcs.r --data $phenofile --pcs $eigenvec --out $newfilepheno $covar2 --npcs $npcs --covar_type $covariates_bin
     pcs=`cat ${newfilepheno}.covar`
     cov_type=`cat ${newfilepheno}.covartype`
     """
}

process format_pheno {
  label 'R'
  input :
       path(phenofile)
       tuple path(bed), path (bim), path(fam)
       val(pheno)
       val(covar)
       val(covar_type)
       val(gxe)
  publishDir "${params.output_dir}/pheno/", mode:'copy'
  output :
       path("${newdatafile}.pheno"), emit :data
       env newpheno, emit : pheno
       val(newcovar), emit : covar
       val(covar_type),emit : covar_type
       val(gxe), emit : gxe
       path("${newdatafile}*")
  script :
      newdatafile=phenofile.baseName+'_format'
      covar2=(covar=="")? "" : " --covar $covar "
      phenores_tr_fct= (params.phenores_tr_fct=="") ? "none"  : "$params.phenores_tr_fct"
      pheno_tr_fct= (params.pheno_tr_fct=="") ? "none"  : "$params.pheno_tr_fct"
      newcovar=(params.pheno_residuals==1) ? "" : " $covar"
      covar_type=(params.pheno_residuals==1) ? "" : " $covar_type"
      gxecmd=(gxe=='') ? "" : "--gxe gxe"  
  """
  format_pheno.r --pheno ${pheno} --data $phenofile $covar2 --out $newdatafile --fam $fam --transform_i ${pheno_tr_fct} --transform_r ${phenores_tr_fct} --residuals ${params.pheno_residuals} $gxecmd
  newpheno=`cat pheno.txt`
  """
}
/*to do check type of covariates*/
workflow wf_prepare_pheno{
   take :
        phenofile
        pheno
        pheno_type
        covar
        covar_type
        gxe
        plink
        pcs
   main :
   if(!pheno)pheno=channel.from(params.pheno)
   if(!covar)covar=channel.from(params.covariates)
   showPhenoDistrib(phenofile,pheno)
   covar_type=covar_type
   if(params.add_pcs>=1){
      add_pcs(pcs.combine(phenofile).combine(channel.from(params.add_pcs)).combine(covar).combine(covar_type))
      new_data=add_pcs.out.data
      newcovar=add_pcs.out.pcs_covar
      data_ch=add_pcs.out.data
      covar_type = add_pcs.out.covar_type
    }else {
    newcovar=covar
    data_ch=phenofile
   }
  format_pheno(data_ch, plink,pheno,newcovar, covar_type, gxe)
  emit :
     data=format_pheno.out.data
     pheno=format_pheno.out.pheno
     covar=format_pheno.out.covar
     covar_type=format_pheno.out.covar_type
     gxe= format_pheno.out.gxe
     plotpheno=showPhenoDistrib.out
}
