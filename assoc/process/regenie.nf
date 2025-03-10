process regenie_step1{
   cpus params.regenie_num_cores
   memory { strmem(params.regenie_mem_req) + 5.GB * (task.attempt -1) }
   input :
     tuple val(pheno),val(covariates),val(covoption_regenie),val(pheno_bin),path(data)
     tuple val(pheno) ,path(bed), path(bim), path(fam), path(rsrel)
   publishDir "${params.output_dir}/assoc/regenie/step1", overwrite:true, mode:'copy', pattern: "*.loco"
   publishDir "${params.output_dir}/assoc/regenie/report/", overwrite:true, mode:'copy', pattern: "*.list"
   publishDir "${params.output_dir}/assoc/regenie/report/", overwrite:true, mode:'copy', pattern: "*.report"
   output :
    tuple val(our_pheno), path("$phef"),path("${out}_pred.list"), path("${out}_1.loco"),path(bed), path(bim), path(fam),  optional :true 
    tuple path("*.report"), path("*.log"), emit : report 
   script :
      phef=pheno+".pheno"
      loco=(params.loco==0) ? "" : " --loocv "
      bfile=bed.baseName
      keeppos=(rsrel.toString()=='00') ? ""   : " --extract $rsrel "
      covoption=(covariates=="") ? "" : " --cov_list ${params.covariates}"
      covoption_regenie=(covariates=="") ? "" : " --covarFile $phef ${covariable_regenie} "
      bsize=(params.regenie_bsize_step1==0) ? " ${params.regenie_bsize} " : " ${params.regenie_bsize_step1}"
      binpheno = (pheno_bin==0) ? "" : " --bt "
      out=phef+"_regenie"
      gxe=(params.regenie_gxe==0) ? "" : " --gxe ${params.gxe} "
      keeppos=(rsrel.toString()=='00') ? ""   : " --extract $rsrel "
      """
      all_covariate.py --data  $data --inp_fam  $fam  $covoption --pheno $pheno --phe_out ${phef}  --form_out 2 --nona  1 $gxe
      ${params.regenie_bin} --step 1   --bed $bfile --phenoFile $phef  --phenoCol ${our_pheno} --bsize $bsize $loco --out  $out --threads ${params.regenie_num_cores} ${params.regenie_otheropt_step1} $covoption_regenie ${binpheno}
      if [ ! -f $out"_1.loco" ]
      then
      touch $out"_1.loco"
      fi
      cp .command.sh "${our_pheno}"_regenie_step1.cmd.report
      cp .command.log "${our_pheno}"_regenie_step1.log.report
      cp .command.err "${our_pheno}"_regenie_step1.err.report
      """
}

 process regenie_step2{
   cpus params.regenie_num_cores
   memory { strmem(params.regenie_mem_req) + 5.GB * (task.attempt -1) }
   input :
    tuple val(pheno), path(data),path(list), path(loco),path(bed), path(bim), path(fam),path(bgen), path(bgensample)
   publishDir "${params.output_dir}/regenie/step2", overwrite:true, mode:'copy', pattern: "*.cmd"
   output :
     tuple val(pheno), path("${out}*${pheno}.regenie")
     path("*.report"),emit : report
   script :
    bfile=bed.baseName
    covoption_regenie= (params.covariates=="") ? "" : " --covarFile $data   ${covariable_regenie} "
    bsize=(params.regenie_bsize_step2==0) ? " ${params.regenie_bsize} " : " ${params.regenie_bsize_step2}"
    genet=(params.bgen=="")? " --bed $bfile " : " --bgen $bgen --sample $bgensample "
    genet=(params.list_bgen=="")? " $genet " : " --bgen $bgen --sample $bgensample "
    loco=(params.loco==0) ? "" : " --loocv "
    out=pheno+"_regenie_assoc"
    out=(params.list_bgen=="")? out: bgen.baseName
    gxe=(params.regenie_gxe==0) ? "" : " --interaction ${params.gxe} "
    """
     ${params.regenie_bin} --step 2  $genet   --phenoFile $data --phenoCol $pheno ${covoption_regenie} --bsize $bsize   --pred $list  $loco   --out $out --threads ${params.regenie_num_cores} ${params.regenie_otheropt_step2} $gxe
      cp .command.sh "${pheno}"_regenie_step2.cmd.report
      cp .command.log "${pheno}"_regenie_step2.log.report
      cp .command.err "${pheno}"_regenie_step2.err.report
  """
 }

 process format_regeniesumstat{
  input :
   tuple val(bfile), val(pheno), path(assoc) 
  publishDir "${params.output_dir}/regenie/", overwrite:true, mode:'copy'
  output :
    tuple val(bfile), val(pheno), path(newassoc) 
  script :
    newassoc=assoc.baseName+"_format.regenie"
  """
  format_sumstat_regenie.py $assoc $newassoc
  """
 }
