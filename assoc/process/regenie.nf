include {merge_sumstat;COFACTORS_TYPE;join2channel} from './all.nf'             
include {indexbgen;bgen_formatsample;getchrobgen} from "../../modules/bgen.nf"  
include { strmem } from "../../modules/fct_groovy.nf"                           

process regenie_step1{
   cpus params.max_cpus
   memory { strmem(params.high_memory) + 5.GB * (task.attempt -1) }
   input :
     tuple val(pheno),val(pheno_bin),val(covariates),val(covoption_regenie),path(data),path(bed), path(bim), path(fam)
   publishDir "${params.output_dir}/assoc/regenie/step1", overwrite:true, mode:'copy', pattern: "*.loco"
   publishDir "${params.output_dir}/assoc/regenie/report/", overwrite:true, mode:'copy', pattern: "*.report"
   output :
    tuple val(pheno), path("$phef"),path("${out}_pred.list"), path("${out}_1.loco"),path(bed), path(bim), path(fam),  emit : step1
    path("*.report"), emit : report 
   script :
      phef=pheno+".pheno"
      loco=(params.loco==0) ? "" : " --loocv "
      bfile=bed.baseName
      covoption=(covariates=="") ? "" : " --cov_list ${params.covariates}"
      covoption_regenie=(covariates=="") ? "" : " --covarFile $phef ${covoption_regenie} "
      bsize=(params.regenie_bsize_step1==0) ? " ${params.regenie_bsize} " : " ${params.regenie_bsize_step1}"
      println(pheno_bin)
      binpheno = (pheno_bin.toInteger()==0) ? "" : " --bt "
      out=phef+"_regenie"
      gxe=(params.gxe=='') ? "" : " --gxe ${params.gxe} "
      """
      all_covariate.py --data  $data --inp_fam  $fam  $covoption --pheno $pheno --phe_out ${phef}  --form_out 2 --nona  1 $gxe
      ${params.regenie_bin} --step 1   --bed $bfile --phenoFile $phef  --phenoCol ${pheno} --bsize $bsize $loco --out  $out --threads ${params.max_cpus} ${params.regenie_otheropt_step1} $covoption_regenie ${binpheno}
      if [ ! -f $out"_1.loco" ]
      then
      touch $out"_1.loco"
      fi
      cp .command.sh "${pheno}"_regenie_step1.cmd.report
      cp .command.log "${pheno}"_regenie_step1.log.report
      cp .command.err "${pheno}"_regenie_step1.err.report
      """
}

 process regenie_step2{
   cpus params.max_cpus
   memory { strmem(params.high_memory) + 5.GB * (task.attempt -1) }
   input :
    tuple val(pheno),path(data),path(list), path(loco),path(bed), path(bim), path(fam),val(covariable_regenie),path(bgen), path(bgensample),val(balise_bgen),val(balise_bgenlist)
   publishDir "${params.output_dir}/regenie/step2", overwrite:true, mode:'copy', pattern: "*.cmd"
   output :
     tuple val(pheno), path("${out}*${pheno}.regenie")
     path("*.report"),emit : report
   script :
    bfile=bed.baseName
    covoption_regenie= (covariable_regenie=="") ? "" : " --covarFile $data   ${covariable_regenie} "
    bsize=(params.regenie_bsize_step2==0) ? " ${params.regenie_bsize} " : " ${params.regenie_bsize_step2}"
    genet=(balise_bgen)? " --bed $bfile " : " --bgen $bgen --sample $bgensample "
    genet=(balise_bgenlist)? " $genet " : " --bgen $bgen --sample $bgensample "
    loco=(loco==0) ? "" : " --loocv "
    out=pheno+"_regenie_assoc"
    out=(balise_bgenlist)? out: bgen.baseName
    gxe=(params.gxe=='') ? "" : " --interaction ${params.gxe} "
    """
     ${params.regenie_bin} --step 2  $genet   --phenoFile $data --phenoCol $pheno ${covoption_regenie} --bsize $bsize   --pred $list  $loco   --out $out --threads ${params.max_cpus} ${params.regenie_otheropt_step2} $gxe
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

workflow regenie {
  take :
   data
   pheno
   pheno_bin
   covariates
   covariates_type
   plink
   plink_rel
   vcf
   vcf_balise
   bgen
   bgenlist
   bgen_sample
   bgen_balise
   bgenlist_balise
   listchro
 main :
 COFACTORS_TYPE(covariates, covariates_type,channel.of(1))
 if(bgen_balise.val || bgenlist_balise.val){
   println("regenie performed using bgen")
   bgen_formatsample(data,bgen_sample)
 }else {
 if(vcf_balise.val){
    println("regenie cannot performed using vcf, to convert in bgen see option")
 }
    println("regenie performed using plink ")
 }
//     tuple val(pheno),val(covariates),val(covoption_regenie),val(pheno_bin),path(data),path(bed), path(bim), path(fam), path(rsrel)
 npheno=params.pheno.split(',').size()                                          
 phenol = pheno.flatMap { list -> list.split(',') }  // Groovy-style lambda for splitting
 pheno_binl = pheno_bin.map { list -> list.split(',') }.flatMap { it.size() == 1 ? it.collect { it } * npheno : it }
 combined = join2channel(phenol,pheno_binl)
 //    tuple val(pheno),val(pheno_bin),val(covariates),val(covoption_regenie),path(data),path(bed), path(bim), path(fam)
 tmp=combined.combine(covariates).combine(COFACTORS_TYPE.out).combine(data).combine(plink_rel)
 regenie_step1(combined.combine(covariates).combine(COFACTORS_TYPE.out).combine(data).combine(plink_rel))
 regenie_step2(regenie_step1.out.step1.combine(COFACTORS_TYPE.out).combine(bgen).combine(bgen_formatsample.out.regenie).combine(bgen_balise).combine(bgenlist_balise))
 
/*
 if(vcf_balise.val){
    println("saige performed using vcf")
    checkidd_saige_vcf(data, pheno,pheno_bin,vcf.first(),plink_rel)
    plkinp=checkidd_saige_vcf.out.plink
   phenoclean=checkidd_saige_saige.out.pheno
 }else if(bgen_balise.val || bgenlist_balise.val){
    println("saige performed using bgen")
   bgen_formatsample(data,bgen_sample)
   checkidd_saige(data, pheno, pheno_bin,bgen_sample, plink_rel)
   plkinp=checkidd_saige.out.plink
   phenoclean=checkidd_saige.out.pheno
 }else{
   checkidd_saige(data, pheno, pheno_bin,channel.fromPath("00"), plink_rel)
   plkinp=checkidd_saige.out.plink
   phenoclean=checkidd_saige.out.pheno
 }*/
}
