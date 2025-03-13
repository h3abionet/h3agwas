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
      tuple path(eigenvec), path(phenofile)
  publishDir "${params.output_dir}/pheno/pca/", mode:'copy'
  output :
    path(newfilepheno), emit :data
    env pcs, emit : pcs_covar
  script :
     plkf=bed.baseName
     head=phenofile.baseName
     newfilepheno=head+"_"+params.add_pcs+".newpheno"
     covar2=(params.covariates=="")? "" : " --covar $params.covariates "

     """
     addpheno_pcs.r --data $phenofile --pcs $head".eigenvec" --out $newfilepheno $covar2
     pcs=`cat ${newfilepheno}.covar`
     """
}

process format_pheno {
  label 'R'
  input :
       path(phenofile)
       tuple path(bed), path (bim), path(fam)
       val(pheno)
       val(covar)
       val(gxe)
  publishDir "${params.output_dir}/pheno/", mode:'copy'
  output :
       path("${newdatafile}.pheno"), emit :data
       env newpheno, emit : pheno
       val(newcovar), emit : covar
       val(gxe), emit : gxe
       path("${newdatafile}*")
  script :
      newdatafile=phenofile.baseName+'_format'
      covar2=(covar=="")? "" : " --covar $covar "
      phenores_tr_fct= (params.phenores_tr_fct=="") ? "none"  : "$params.phenores_tr_fct"
      pheno_tr_fct= (params.pheno_tr_fct=="") ? "none"  : "$params.pheno_tr_fct"
      newcovar=(params.pheno_residuals==1) ? "" : " $covar"
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
        covar
        gxe
        plink
        pcs
   main :
   if(!pheno)pheno=channel.from(params.pheno)
   if(!covar)covar=channel.from(params.covariates)
   println(params.covariates)
   showPhenoDistrib(phenofile,pheno)
   if(params.add_pcs>=1){
      add_pcs(pcs.combine(phenofile))
      new_data=add_pcs.out.data
      newcovar=add_pcs.out.pcs_covar
      //if(params.covar!='')newcovar=params.covar+","+newcovar
      data_ch=add_pcs.out.data
    }else {
    newcovar=params.covariates
    data_ch=phenofile
   }
  format_pheno(data_ch, plink,pheno,newcovar, gxe)
  emit :
     data=format_pheno.out.data
     pheno=format_pheno.out.pheno
     covar=format_pheno.out.covar
     gxe= format_pheno.out.gxe
     plotpheno=showPhenoDistrib.out
}
