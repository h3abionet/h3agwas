#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Shaun Aron
 *      Rob Clucas
 *      Eugene de Beste
 *      Lerato Magosi
 *      Brandenburg Jean-Tristan
 *
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS and simulation data
 *
 */


//---- General definitions --------------------------------------------------//

import java.nio.file.Paths


def helps = [ 'help' : 'help' ]

allowed_params = ["data","covariates","fastlmm_multi","max_plink_cores","gemma_num_cores","other_mem_req","plink_mem_req","input_dir","input_pat","output","output_dir","num_cores","mem_req","gemma","linear","logistic","chi2","fisher", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "sharedStorageMount", "shared-storage-mount",  "big_time","thin"]

param_bolt=["bolt_ld_scores_col","bolt_ld_score_file","boltlmm", "bolt_covariates_type",  "bolt_use_missing_cov"]
allowed_params+=param_bolt
param_phenosim=["ph_cov_norm","ph_mem_req","ph_cov_range","ph_intercept","ph_nb_qtl","ph_list_qtl","ph_nb_sim", "ph_quant_trait", "ph_qual_dom", "ph_maf_r", "ph_alpha_lim", "ph_windows_size", "ph_normalise"]
allowed_params+=param_phenosim


params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}

def params_help = new LinkedHashMap(helps)

params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output_testing = "cleaned"
params.thin       = ""
params.covariates = ""
params.chrom      = ""
outfname = params.output_testing
params.gemma_num_cores = 6
params.big_time="100H"

supported_tests = ["chi2","fisher","model","cmh","linear","logistic","boltlmm", "fastlmm", "gemma"]


params.chi2     = 0
params.fisher   = 0
params.cmh     =  0
params.model   =  0
params.linear   = 0
params.logistic = 0
params.gemma = 0
params.mem_req = "6GB"
params.gemma_relopt = 1
params.gemma_lmmopt = 4

/*JT Append initialisation variable*/
params.bolt_covariates_type = ""
params.bolt_ld_score_file= ""
params.bolt_ld_scores_col=""
params.boltlmm = 0
params.bolt_use_missing_cov=0
params.num_cores=1

params.input_pat  = 'raw-GWA-data'

params.sexinfo_available = "false"
params.data=""

/*param for phenosim*/
/*Number simulation*/
params.ph_nb_sim=5
/*quantitative traits => 1
qualitative traits 0*/
params.ph_quant_trait=1
/*Qualitative */
/*Nb qtl*/
params.ph_nb_qtl=2
/*qtl : number be as ph_nb_qtl, separate by ","*/
params.ph_list_qtl=""
/*ph_qual_dom Quantitative : adititve : 0 dominant 1 */
params.ph_qual_dom=0
/*freq for each snps*/
params.ph_maf_r="0.05,1.0"
params.ph_alpha_lim="0.05,0.000001"
params.ph_windows_size="1000000bp"
/*balise for normalisation */
/* used variable pheno, covar to normalise with specific */
params.ph_normalise = 0
/*need to be used if params normalised */
/*cov1=x,cov2=y*/
params.ph_eq_norm=""
params.ph_cov_norm=""
params.ph_intercept=""
params.ph_cov_range=""
params.ph_mem_req="20GB"

params.help = false
if (params.help) {
    params.each {
    entry ->
      print "Parameter: <$entry.key>    \t Default: $entry.value"
      if (entry.key == 'help')
          println ""
      else {
        help = params_help.get(entry.key)
        if (help)
          print "\n    $help"
        println ""
      }
  }
  System.exit(-1)
}

def getConfig = {
  all_files = workflow.configFiles.unique()
  text = ""
  all_files.each { fname ->
      base = fname.baseName
      curr = "\n\n*-subsection{*-protect*-url{$base}}@.@@.@*-footnotesize@.@*-begin{verbatim}"
      file(fname).eachLine { String line ->
        if (line.contains("secretKey")) { line = "secretKey='*******'" }
        if (line.contains("accessKey")) { line = "accessKey='*******'" }
        curr = curr + "@.@"+line
      }
      curr = curr +"@.@*-end{verbatim}\n"
      text = text+curr
  }
  return text
}


// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n------\nError in your config\nFile $fn does not exist\n\n---\n")
}


if (params.thin)
   thin = "--thin ${params.thin}"
else
   thin = ""

if (params.chrom)
   chrom = "--chr ${params.chrom}"
else
   chrom = ""

if(params.ph_normalise==1){
if(params.data=="")error("ph_normalise must be have a data file")
}
/*Initialisation of bed data, check files exist*/
raw_src_ch= Channel.create()


bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


/*if need a chromosome or sub sample of data*/
bed_all_file_ms=Channel.create()
if (thin+chrom) {
  process thin {
    time params.big_time
    input:
      set file(bed), file(bim), file(fam) from raw_src_ch
    output:
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into (bed_all_file_msI)
    script:
       base = bed.baseName
       out  = base+"_t"
       "plink --keep-allele-order --bfile $base $thin $chrom --make-bed --out $out"
  }
  }else{
        bed_all_file_msI=Channel.create()
        raw_src_ch.into(bed_all_file_msI)
 }
if(params.ph_normalise){
  liste_cov=""
  if(params.ph_cov_norm!=""){
   liste_cov=params.ph_cov_norm.replaceAll(/=[0-9.-]+,/,",").replaceAll(/=[0-9.-]+$/,"")+","
  }
  if(params.covariates!=""){
    liste_cov+=params.covariates
  }
  liste_cov=liste_cov.replaceAll(/,$/,"")
  data_ch = Channel.fromPath(params.data)
  if(liste_cov=="")liste_cov="FID"
  process GetDataNoMissing{
     time params.big_time
     input :
       set file(bed), file(bim), file(fam) from bed_all_file_msI
       file(data) from data_ch
     output :
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into (bed_all_file_ms)
     script :
       base = bed.baseName
       out = base+"nona"
       data_nomissing=".temp"
       list_ind_nomissing="listeInd"
       """
       list_ind_nomissing.py --data $data --inp_fam $fam --cov_list $liste_cov --pheno FID --dataout  $data_nomissing --lindout $list_ind_nomissing 
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $out
       """
}
}else{
  bed_all_file_ms=Channel.create()
  bed_all_file_msI.into(bed_all_file_ms)
}
/*transform data in ms for phenosim*/
process ChangeMsFormat{
   time params.big_time
   cpus params.num_cores
   memory params.mem_req
   input :
     set file(bed), file(bim), file(fam) from bed_all_file_ms
   output :
     set file(data_ms), file(data_ms_ped), file(data_ms_bed), file(data_ms_bim), file(data_ms_fam) into phenosim_data 
     set file(data_ms_bed), file(data_ms_bim), file(data_ms_fam) into data_gemma_rel
   script :
     base=bed.baseName
     base_ped=base+"_ped"
     data_ms=base+".ms"
     data_ms_ped=data_ms+".ped"
     data_ms_map=data_ms+".map"
     data_ms_fam=data_ms+".fam"
     data_ms_bim=data_ms+".bim"
     data_ms_bed=data_ms+".bed"
     """
     plink --keep-allele-order --bfile $base --threads ${params.num_cores} --recode tab --out $base_ped
     convert_plink_ms.py $base_ped".ped" $bim $data_ms $data_ms_ped
     cp $base_ped".map"  $data_ms_map
     cp $base".fam"  $data_ms_fam
     cp $base".bim"  $data_ms_bim
     plink --keep-allele-order --file $data_ms --make-bed  --out $data_ms --threads ${params.num_cores}
     """
}

phenosim_data_all=Channel.from(1..params.ph_nb_sim).combine(phenosim_data)
if(params.ph_quant_trait==1){
if(!params.ph_list_qtl)LQTL=([0.05]*params.ph_nb_qtl).join(",")
else LQTL=params.ph_list_qtl
ph_quant_trait_param="-n ${params.ph_nb_qtl} -v $LQTL"
}
process SimulPheno{
   memory params.ph_mem_req
   time params.big_time
   input :
     set sim, file(ms), file(ped),file(bed), file(bim), file(fam) from phenosim_data_all
   publishDir "${params.output_dir}/simul/", pattern: "$ent_out_phen*", overwrite:true, mode:'copy'
   output :
     set sim, file(file_causal), file(file_pheno), file(bed), file(bim), file(fam) into raw_pheno
   script :
     base=bed.baseName
     ent_out_phen=base+"."+sim+".pheno"
     file_pheno=ent_out_phen+".pheno"
     file_causal=ent_out_phen+".causal"
     """
     phenosim.py -d 1 -f $ms -i M -o N -q ${params.ph_quant_trait} $ph_quant_trait_param --maf_c 0 --outfile $ent_out_phen
     echo -e "FID\\tIID\\tPhenoS" > $file_pheno
     paste $fam $ent_out_phen"0.pheno"|awk '{print \$1"\\t"\$2"\\t"\$NF}' >> $file_pheno
     ListePos=`awk 'BEGIN{AA=""}{AA=AA\$2","}END{print AA}' ${ent_out_phen}0.causal`
     get_lines_bynum.py --file $bim --lines \$ListePos --out ${file_causal}"tmp"
     paste ${file_causal}"tmp" ${ent_out_phen}0.causal > $file_causal
     """
}
if(params.ph_normalise && params.data){
   normal_cov_param=""
   if(params.ph_cov_norm){
      normal_cov_param+="--cov_info="+params.ph_cov_norm
      normal_cov_param+=" --rangenorm="+params.ph_cov_range
      if(params.ph_intercept)normal_cov_param+=" --intercept="+params.ph_intercept
   }
   raw_pheno_data= Channel.fromPath(params.data).combine(raw_pheno)
   process NormaliseDataph_normalise{
     cpus params.num_cores
     memory params.mem_req
     time params.big_time
     input :
       set file(data), sim, file(file_causal), file(file_pheno), file(bed), file(bim), file(fam) from raw_pheno_data
     output :
       set sim, file(file_causal), file(file_pheno_new), file(bed), file(bim), file(fam) into (sim_data_gemma,sim_data_bolt, sim_data_plink)
     script :
       file_pheno_new=file_pheno+".norm"
       """
       ph_normalise_variable.py --data $data $normal_cov_param --data_sim $file_pheno --out $file_pheno_new --na_out "NA"
       """
   }
}else{
  sim_data_gemma=Channel.create()
  sim_data_bolt=Channel.create()
  sim_data_plink=Channel.create()
  raw_pheno.into(sim_data_gemma,sim_data_bolt, sim_data_plink)
}
if(params.gemma==1){
  process getGemmaRel {
    cpus params.num_cores
    memory params.mem_req
    time params.big_time
    input:
       file plinks from data_gemma_rel
    output:
       file("output/${base}.*XX.txt") into rel_mat_ch
    script:
       base = plinks[0].baseName
       famfile=base+".fam"

       """
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       cat $famfile |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
       gemma -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3 
       """
  }
  sim_data_gemma2=sim_data_gemma.combine(rel_mat_ch)
  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""
  process doGemma{
    cpus params.num_cores
    memory params.mem_req
    time params.big_time
    input:
      set sim, file(file_causal), file(file_pheno), file(bed), file(bim), file(fam), file(rel) from sim_data_gemma2
    publishDir "${params.output_dir}/gemma/simul/", pattern: "output/${out}.assoc.txt", overwrite:true, mode:'copy'
    output :
      set sim, file(file_causal), file("output/${out}.assoc.txt"), file(bim) into res_gem
    script :
       base = bed.baseName
       gemma_covariate=base+"."+sim+".cov"
       gemma_pheno=base+"."+sim+".phe"
       out=base+"."+sim+".gem"
       covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
       """
       all_covariate.py --data  $file_pheno --inp_fam  $fam ${covariate_option} --cov_out $gemma_covariate \
                          --pheno PhenoS --phe_out $gemma_pheno --form_out 1
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       gemma -bfile $base ${covar_opt_gemma}  -k $rel -lmm 1  -n 1 -p $gemma_pheno  -o $out 
       """
   }
   alpha_lim=params.ph_alpha_lim
   process doGemmaStat{
    input:
      set sim, file(file_causal), file(stat), file(bim) from res_gem
    output :
     file(out) into res_stat_gem
    script :
      base=bim.baseName
      out=base+".gem."+sim+".res.stat"
      """
      compute_stat_phenosim.py --stat $stat --bim  $bim --header_pval p_wald --header_chro chr --header_pos ps --windows_size $params.ph_windows_size --alpha_lim $params.ph_alpha_lim --out $out --pos_simul $file_causal
      """
  }
  res_stat_gemmac=res_stat_gem.collect()
  process doMergeStatGemma{
    input :
     val(liste_file) from res_stat_gemmac
    publishDir "${params.output_dir}/gemma", overwrite:true, mode:'copy'
    output :
     file(output_gemma)
    script :
      output_gemma="res_gemma.stat"
      listefiles=liste_file.join(" ")
      file = liste_file[0]
      """
      head -1 $file > $output_gemma
      cat $listefiles|grep -v "nsig" >> $output_gemma
      """
  }


}
   /*JT Fonction to transforme argument for cofactor in gemma
   @Input
   args: cofactor args separate by a comma
   infoargs: type of cofactor separate by a comma : 0 for qualitative, 1 for quantitative
   output : cofactor for boltlmm was formating to take account qualitative and quantitative
   */
   def boltlmmCofact(args,infoargs) {
      //Method code
      splargs=args.split(",")
      splinfoargs=infoargs.split(",")
      if(splargs.size() != splinfoargs.size()){
         System.err.println("args and args type for Boltlmm was not same size : "+args+" "+infoargs)
         System.exit(-11)
      }
      CofactStr=""
      for (i = 0; i <splargs.size(); i++) {
          /*0 : for quantitatif */
          /* 1 for qualitatif*/
          if     (splinfoargs[i]=='1')  CofactStr +=" --qCovarCol="+splargs[i]
          else if(splinfoargs[i]=='0')  CofactStr+=" --covarCol="+splargs[i]
          else{
             System.err.println("type args for "+splargs[i]+" doesn't know "+ splinfoargs[i]+"\n 1 for quantitatif arguments\n 0 for qualitatif arguments")
             System.exit(-10)
          }
      }
      return(CofactStr)
   }


if(params.boltlmm==1){
   if(params.bolt_ld_score_file){
      ld_score_cmd="--LDscoresFile=$params.bolt_ld_score_file"
      if(params.bolt_ld_scores_col) ld_score_cmd = ld_score_cmd + " --LDscoresCol="+params.bolt_ld_scores_col
   } else
      ld_score_cmd = "--LDscoresUseChip"

   if (params.covariates)
      bolt_covparam= boltlmmCofact(params.covariates,params.bolt_covariates_type)
   else
      bolt_covparam= ""

  type_lmm="--lmm"
  
  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""

  process doBoltlmmm{
    cpus params.num_cores
    memory params.mem_req
    time params.big_time
    input:
      set sim, file(file_causal), file(file_pheno), file(bed), file(bim), file(fam) from sim_data_bolt
    publishDir "${params.output_dir}/boltlmm/simul/", pattern: "${out}", overwrite:true, mode:'copy'
    output :
      set sim, file(file_causal), file(out), file(bim) into res_bolt
    script :
       base = bed.baseName
       bolt_covariate=base+"."+sim+".cov"
       bolt_pheno=base+"."+sim+".phe"
       out=base+"."+sim+".bolt"
       covar_file_bolt =  (params.covariates) ?  " --covarFile ${bolt_pheno} " : ""
       """
       shuf -n 950000 $bim | awk '{print \$2}' > .sample.snp
       all_covariate.py --data  $file_pheno --inp_fam  $fam ${covariate_option} --cov_out $bolt_covariate \
                          --pheno PhenoS --phe_out $bolt_pheno --form_out 2
       bolt $type_lmm --bfile=$base  --phenoFile=$bolt_pheno --phenoCol=PhenoS --numThreads=$params.num_cores  $covar_file_bolt --statsFile=$out\
           $ld_score_cmd  --lmmForceNonInf  --modelSnps=.sample.snp $bolt_covparam
       """
   }
   alpha_lim=params.ph_alpha_lim
   process doBoltStat{

    input:
      set sim, file(file_causal), file(stat), file(bim) from res_bolt
    output :
      file(out) into res_stat_bolt
    script :
      base=bim.baseName
      out=base+"."+sim+".bolt.res.stat"
      """
      compute_stat_phenosim.py --stat $stat --bim  $bim --header_pval P_BOLT_LMM_INF --header_chro CHR --header_pos BP --windows_size $params.ph_windows_size --alpha_lim $params.ph_alpha_lim --out $out --pos_simul $file_causal
      """

  }
   res_stat_boltc=res_stat_bolt.collect()
  process doMergeStatBolt{
    input :
      val(liste_file) from res_stat_boltc
    publishDir "${params.output_dir}/boltlmm", overwrite:true, mode:'copy'
    output :
       file(output_bolt) 
    script :
      output_bolt="res_boltlmm.stat"
      listefiles=liste_file.join(" ")
      file = liste_file[0]
      """
      head -1 $file > $output_bolt
      cat $listefiles|grep -v "nsig" >> $output_bolt
      """ 
  }
}







