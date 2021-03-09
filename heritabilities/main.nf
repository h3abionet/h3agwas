#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description : pipeline to estimate heritabilies and co heritabilites with summary stat, genetics data.... with gemma, bolt, ldsc, gcta
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths


def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin", "gc10"]

// define param for
annotation_param=[ "file_gwas", "Nind"]
allowed_params+=annotation_param
h2_param=[ "bolt_h2", "gcta_h2", "gcta_h2_imp","bolt_h2_multi"]
allowed_params+=h2_param
h2_bolt=["bolt_ld_scores_col", "bolt_ld_score_file","boltlmm", "bolt_covariates_type",  "bolt_use_missing_cov", "bolt_num_cores", "bolt_mem_req", "exclude_snps", "bolt_impute2filelist", "bolt_impute2fidiid", "bolt_otheropt"]
h2_gemma=["gemma_h2_pval","gemma_num_cores", "gemma_mem_req", "gemma_h2", "gemma_h2_pval", "gemma_h2_typeest"]
h2_ldsc=["ldsc_h2", "ldsc_h2_multi", "ldsc_mem_req","ldsc_h2opt", "genetic_map_file","dir_ref_ld_chr", "ldsc_bin","munge_sumstats_bin"]
h2_gcta=["gcta_h2_multi", "grm_cutoff","multigrm_opt", "gcta_num_cores", "gcta_bin", "gcta_mem_req"]
allowed_params+=h2_gemma
allowed_params+=h2_bolt
allowed_params+=h2_gcta
allowed_params+=h2_ldsc


allowed_params_head = ["head_pval", "head_freq", "head_bp", "head_chr", "head_rs", "head_beta", "head_se", "head_A1", "head_A2", "head_n"]
allowed_params+=allowed_params_head

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
params.covariates = ""
params.Nind=0
params.pheno=""
params.munge_sumstats_bin="munge_sumstats.py"
params.ldsc_bin="ldsc.py"
outfname = params.output_testing
params.file_rs_buildrelat = ""
params.genetic_map_file=""

params.file_gwas=""
params.head_pval = "P_BOLT_LMM"
params.head_freq = "A1FREQ"
params.head_n = ""
params.head_bp = "BP"
params.head_chr = "CHR"
params.head_rs = "SNP"
params.head_beta="BETA"
params.head_se="SE"
params.head_A1="ALLELE1"
params.head_A2="ALLELE0"

//gcta grm option 
params.grm_cutoff=0.025
params.multigrm_opt=""
params.gcta_num_cores = 10

// bolt option default
params.bolt_h2=0
params.bolt_covariates_type = ""
params.bolt_ld_score_file= ""
params.bolt_ld_scores_col=""
params.bolt_num_cores=8
params.bolt_mem_req="6GB"
params.bolt_use_missing_cov=0
params.exclude_snps=""
params.bolt_impute2filelist=""
params.bolt_impute2fidiid=""
params.bolt_otheropt=""
params.bolt_h2_multi=0
params.gcta_grmfile=""



params.cut_maf=0.01
//pos with rsid pos chro
params.list_snp=""
//## info if need
params.cut_info=0.6
params.plink_mem_req="6GB"
params.data=""


params.big_time='1000h'
//ldsc
params.ldsc_h2=0
params.ldsc_h2_multi=0
params.ldsc_mem_req="6GB"
params.ldsc_h2opt=""
params.dir_ref_ld_chr=""
//params for gcta
params.gcta_bin="gcta64"
params.gcta_h2=0
params.gcta_h2_multi=0
params.gcta_mem_req="15GB"
params.gcta_mem_reqmgrm = "40GB"
params.gcta_h2_imp=0
params.gcta_h2_ldscore = 200
params.gcta_h2_grmfile =""
params.gcta_reml_alg=0
params.gcta_opt_multigrm_cor=""
params.gcta_opt_grm_cor=""

//params gemma
params.gemma_bin="gemma"
params.gemma_h2=0
params.gemma_h2_pval = 0
/*1 ou 2*/
params.gemma_h2_typeest="1"
params.gemma_mat_rel=""
params.gemma_num_cores=6
params.gemma_mem_req="10GB"
params.gemma_relopt = 1



params.output="heritabilies"

gcta_mem_req=params.gcta_mem_req
gcta_num_cores = params.gcta_num_cores+1
plink_mem_req = params.plink_mem_req
max_plink_cores = params.max_plink_cores
ldsc_mem_req=params.ldsc_mem_req
params.help = false

if(params.gcta_h2_grmfile=="" && params.gcta_grmfile!=""){
params.gcta_h2_grmfile=params.gcta_grmfile
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

println "\nTesting data : ${params.data}\n"
println "Testing for gwas file : ${params.file_gwas}\n"

//////
///
///
///  LDSC
///
//////

if(params.ldsc_h2==1){
 println "ldsc_h2"
 if(params.file_gwas==""){
  error("\n\n------\nbCan't do ldsc without file_gwas\n\n---\n")
 }
 gwas_file_ldsc=Channel.from(params.file_gwas.split(",")).flatMap{it->file(it ,checkIfExists:true)}
 if(params.head_n=="" && params.Nind==0){
       error("\n\n------\nheader_n or Nind paramaters must be allowed\n\n---\n")
 }
 if(params.dir_ref_ld_chr==""){
       error("\n\n------\ndir_ref_ld_chr must be allowed with full path\n\n---\n")

 }
 process DoLDSC{
   memory ldsc_mem_req
   input :
      file(gwas) from gwas_file_ldsc
   publishDir "${params.output_dir}/ldsc", overwrite:true, mode:'copy'
   output :
     file("$out"+".log") into report_ldsc
   script :
     NInfo=params.head_n=="" ? " --N ${params.Nind} " : "--N-col ${params.head_n} " 
     out=gwas.baseName.replace('_','-')+"_ldsc"
     """
     ${params.munge_sumstats_bin} --sumstats $gwas $NInfo --out $out"_mg" --snp ${params.head_rs} --p ${params.head_pval} \
     --frq ${params.head_freq} --info-min ${params.cut_info} --maf-min ${params.cut_maf} --a1 ${params.head_A1} --a2 ${params.head_A2}  --no-alleles
     ${params.ldsc_bin} --h2 $out"_mg.sumstats.gz" --ref-ld-chr ${params.dir_ref_ld_chr} --w-ld-chr ${params.dir_ref_ld_chr} --out $out ${params.ldsc_h2opt}
     format_correlation.r $out".log" $out".stat" ldsc $gwas None
     """ 
}
}else{
 report_ldsc=Channel.empty()
}
if(params.ldsc_h2_multi==1){
 println "ldsc_h2_multi"
 if(params.file_gwas==""){
  error("---------------\n ldsc_h2_multi=1 and file_gwas is null\n\n---------------")
 }
 gwas_file_1 = Channel.fromPath(params.file_gwas.split(",")[0], checkIfExists:true)
 if(params.list_snp==""){
  process GetSnpList{
   time params.big_time
   input :
      file(gwas) from gwas_file_1
   output :
     file("$out") into list_file_pos
   script :
     out="info_poschro.list"
     """
     ma_extract_rsid.py --input_file $gwas --out_file $out --info_file rsID:${params.head_rs},A1:${params.head_A1},A2:${params.head_A2} --ldsc
     """
  }
 }else{
  list_file_pos=Channel.fromPath(params.list_snp,checkIfExists:true) 
 }

  list_gwas_multi=Channel.from(params.file_gwas.split(",")).flatMap{it->file(it, checkIfExists:true)}.combine(list_file_pos)
  process FormatLDSC{
   time params.big_time
   memory params.ldsc_mem_req
   memory ldsc_mem_req
   input :
    set file(gwas), file(infopos) from list_gwas_multi
   output :
    file("$out"+".sumstats.gz") into (list_gwas_formatldsc, list_gwas_formatldsc2)
   script :
    out=gwas.baseName.replace('_','-')
    NInfo=params.head_n=="" ? " --N ${params.Nind} " : "--N-col ${params.head_n} " 
    """
    ${params.munge_sumstats_bin} --sumstats $gwas $NInfo --out $out --snp ${params.head_rs} --p ${params.head_pval} \
    --frq ${params.head_freq} --info-min ${params.cut_info} --maf-min ${params.cut_maf} --a1 ${params.head_A1} --a2 ${params.head_A2}  --merge-alleles $infopos 
    """
  }

 list_gwas_formatldsc3=list_gwas_formatldsc.collect()
 process DoCorrLDSC{
  time params.big_time
  memory ldsc_mem_req
  input :
   file(listfilegwas) from list_gwas_formatldsc3
  publishDir "${params.output_dir}/ldsc", overwrite:true, mode:'copy'
   each pos from 1..params.file_gwas.split(",").size()
  output :
    file("$out"+".log")
  script :
    filegwas=listfilegwas[pos-1]
    listfilegwas.remove((pos-1))
    listfil=filegwas+","+listfilegwas.join(',')
    out = filegwas.baseName.replace('_','-')+"_ldsc_mc"
    println listfil
    println filegwas
    """
    ${params.ldsc_bin} --rg $listfil \
    --ref-ld-chr ${params.dir_ref_ld_chr}\
    --w-ld-chr  ${params.dir_ref_ld_chr} \
    --out  $out
    """ 
  }
}


/*
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
          else if(splinfoargs[i]=='0')  CofactStr +=" --covarCol="+splargs[i]
          else{
             System.err.println("type args for "+splargs[i]+" doesn't know "+ splinfoargs[i]+"\n 1 for quantitatif arguments\n 0 for qualitatif arguments")
             System.exit(-10)
          }
      }
      return(CofactStr)
   }
     def CountLinesFile(File){
     BufferedReader reader = new BufferedReader(new FileReader(File));
     int lines = 0;
     while (reader.readLine() != null) lines++;
     reader.close();
     return(lines)
  }
 
//////
///
///
///  Bolt
///
//////

if(params.bolt_h2){
    println "bolt_h2"
    bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
    bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
    fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
    boltlmm_assoc_ch= Channel.create()
    Channel.from(file(bed),file(bim),file(fam)).buffer(size:3)
        .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
        .set { boltlmm_assoc_ch }

    plink_ch_bolt = Channel.create()
    plink_ch_bolt_multi = Channel.create()
    bim_ch_bolt_snpchoice = Channel.create()
    fam_ch_bolt = Channel.create()
    bim_ch_bolt = Channel.create()
    boltlmm_assoc_ch.separate (plink_ch_bolt, plink_ch_bolt_multi,fam_ch_bolt, bim_ch_bolt, bim_ch_bolt_snpchoice) { a -> [ a, a,a[2], a[1],a[1]] }
    data_ch_bolt = Channel.fromPath(params.data)
    if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
    else
     covariate_option = ""
    process  getBoltPhenosCovar {
      input:
        file(covariates) from data_ch_bolt
        file(fam) from fam_ch_bolt
      output:
        file(phef) into (newdata_ch_bolt,newdata_ch_bolt_multi)
        stdout into (pheno_cols_ch_bolt,pheno_cols_ch_bolt_multi)
      script:
        base = fam.baseName
        phef = "${base}_bolt_n.phe"
        """
        all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
                          --pheno ${params.pheno} --phe_out ${phef} --form_out 2
        """
   }
   ind_pheno_cols_ch_bolt = Channel.create()
   check_bolt = Channel.create()
   pheno_cols_ch_bolt.flatMap { list_str -> list_str.split() }.tap ( check_bolt) .set { ind_pheno_cols_ch_bolt}


   if (params.covariates)
      cov_bolt = boltlmmCofact(params.covariates,params.bolt_covariates_type)
   else
      cov_bolt= ""

   missing_cov=""
   if(params.bolt_use_missing_cov==1)
     missing_cov=" --covarUseMissingIndic "

   pval_head = "P_BOLT_LMM"

  if(params.exclude_snps){
     println "snp exclude of files "+params.exclude_snps
       rs_ch_exclude_bolt=Channel.fromPath(params.exclude_snps)
       rs_ch_exclude_bolt_multi=Channel.fromPath(params.exclude_snps)
  }else{
     println "no snp exclude"
     rs_ch_exclude_bolt=file('NO_FILE')
     rs_ch_exclude_bolt_multi=file('NO_FILE')
  }
  if(params.file_rs_buildrelat!=""){
      println "snp used for rs :"+params.file_rs_buildrelat
      filers_matrel=Channel.fromPath(params.file_rs_buildrelat)
      BoltNbMaxSnps=CountLinesFile(params.file_rs_buildrelat)
      filers_matrel_multi=Channel.fromPath(params.file_rs_buildrelat)
  }else{
      println "CHoice of random Snps to build reladness"
      BoltNbMaxSnps=1000000
      process buildBoltFileSnpRel{
         memory params.bolt_mem_req
         time   params.big_time
         input:
           file(plinksbim) from bim_ch_bolt_snpchoice
         output :
           file(output) into (filers_matrel,filers_matrel_multi)
         script :
           output=plinksbim.baseName.replace('_','-')+".rs.choice"
           """
           shuf -n 950000 $plinksbim | awk '{print \$2}' > $output
           """
      }

  }


  if(params.bolt_ld_score_file!=""){
     println "bolt : ld files used : "+params.bolt_ld_score_file
     Bolt_ld_score= Channel.fromPath(params.bolt_ld_score_file)
     Bolt_ld_score_multi= Channel.fromPath(params.bolt_ld_score_file)
  }else{
     println "no ld files used for bolt "
     Bolt_ld_score = file('NO_FILE3')
     Bolt_ld_score_multi = file('NO_FILE3')
  }
//genetic_map_file
  if(params.genetic_map_file!=""){
     println "genetic map used : "+params.genetic_map_file
     Bolt_genetic_map= Channel.fromPath(params.genetic_map_file)
     Bolt_genetic_map_multi= Channel.fromPath(params.genetic_map_file)
  }else{
     println "no genetic maps used "
     Bolt_genetic_map = file('NO_FILE4')
     Bolt_genetic_map_multi = file('NO_FILE4')
  }


  process doh2Bolt{
    cpus params.bolt_num_cores
    memory params.bolt_mem_req
    time   params.big_time
    input:
      set file(plinksbed), file(plinksbim), file(plinksfam) from plink_ch_bolt
      file(phef) from newdata_ch_bolt
      file(rs_exclude) from rs_ch_exclude_bolt
      file(SnpChoiceMod) from filers_matrel
      file(bolt_ld_score) from Bolt_ld_score
      file(bolt_genetic_map) from Bolt_genetic_map
    publishDir "${params.output_dir}/boltlmm", overwrite:true, mode:'copy'
    each this_pheno from ind_pheno_cols_ch_bolt
    output:
      file(outReml)
      file("${outReml}_bolt.stat") into report_bolt
    script:
      base = plinksbed.baseName
      our_pheno = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"").replaceAll(/^[0-9]+@@@/,"")
      our_pheno2 = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"").replaceAll(/^[0-9]+@@@/,"")
      our_pheno3 = this_pheno.replaceAll(/\/np.\w+/,"").replaceAll(/-$/,"").replaceAll(/^[0-9]+@@@/,"")
      outReml = "$our_pheno2"+".reml"
      covar_file_bolt =  (params.covariates) ?  " --covarFile ${phef} " : ""
      model_snp  = "--modelSnps=$SnpChoiceMod --maxModelSnps=$BoltNbMaxSnps "
      ld_score_cmd = (params.bolt_ld_score_file!="") ? "--LDscoresFile=$bolt_ld_score" :" --LDscoresUseChip "
      ld_score_cmd = (params.bolt_ld_score_file!="" & params.bolt_ld_scores_col!="") ? "$ld_score_cmd --LDscoresCol=${params.bolt_ld_scores_col}" :" $ld_score_cmd "
      exclude_snp = (params.exclude_snps!="") ? " --exclude $rs_exclude " : ""
      geneticmap = (params.genetic_map_file!="") ?  " --geneticMapFile=$bolt_genetic_map " : ""
      """
      bolt.py bolt  --reml  --bfile=$base  --phenoFile=${phef} --phenoCol=${our_pheno3} --numThreads=$params.bolt_num_cores $cov_bolt $covar_file_bolt $missing_cov $model_snp $geneticmap $exclude_snp $ld_score_cmd ${params.bolt_otheropt} --out_bolt2 $outReml
     format_correlation.r $outReml $outReml"_bolt.stat" bolt ${our_pheno} None
      """
  }

if (params.bolt_h2_multi==1){
  println "bolt multi done"
  process doh2BoltiMulti{
    cpus params.bolt_num_cores
    memory params.bolt_mem_req
    time   params.big_time
    input:
      set file(plinksbed), file(plinksbim), file(plinksfam) from plink_ch_bolt_multi
      file(phef) from newdata_ch_bolt_multi
      file(rs_exclude) from rs_ch_exclude_bolt_multi
      file(SnpChoiceMod) from filers_matrel_multi
      file(bolt_ld_score) from Bolt_ld_score_multi
      file(bolt_genetic_map) from Bolt_genetic_map_multi
      val phenoinfo from pheno_cols_ch_bolt_multi
    publishDir "${params.output_dir}/boltlmm_multi", overwrite:true, mode:'copy'
    output:
      file(outReml_multi)
   script :
      phenonew="--phenoCol="+params.pheno.split(',').join(" --phenoCol=")
      base = plinksbed.baseName
      outReml_multi = "$base-all"+".reml"
      covar_file_bolt =  (params.covariates) ?  " --covarFile ${phef} " : ""
      model_snp  = "--modelSnps=$SnpChoiceMod --maxModelSnps=$BoltNbMaxSnps "
      ld_score_cmd = (params.bolt_ld_score_file!="") ? "--LDscoresFile=$bolt_ld_score" :" --LDscoresUseChip "
      ld_score_cmd = (params.bolt_ld_score_file!="" & params.bolt_ld_scores_col!="") ? "$ld_score_cmd --LDscoresCol=${params.bolt_ld_scores_col}" :" $ld_score_cmd "
      exclude_snp = (params.exclude_snps!="") ? " --exclude $rs_exclude " : ""
      geneticmap = (params.genetic_map_file!="") ?  " --geneticMapFile=$bolt_genetic_map " : ""
      """
      bolt.py bolt  --reml  --bfile=$base  --phenoFile=${phef} $phenonew --numThreads=$params.bolt_num_cores $cov_bolt $covar_file_bolt $missing_cov $model_snp $geneticmap $exclude_snp $ld_score_cmd ${params.bolt_otheropt} --out_bolt2 ${outReml_multi}
      """
  }
}
}else{
report_bolt=Channel.empty()
}

//////
///
///
/// GCTA GREML
///
//////

if(params.gcta_h2==1){
   gctabed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
   gctabim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
   gctafam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
   h2gcta_assoc_ch= Channel.create()
   Channel
    .from(file(gctabed),file(gctabim),file(gctafam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { h2gcta_assoc_ch}

  fam_ch_gcta = Channel.create()
  bim_ch_gcta = Channel.create()
  plink_ch_gcta=Channel.create()
  plink_ch_gcta_multigrm=Channel.create()
  plink_ch_gcta_grm=Channel.create()
  fam_ch_gcta_mult = Channel.create()
  plink_ch_gcta_multigrm2 = Channel.create()
  h2gcta_assoc_ch.separate (plink_ch_gcta,plink_ch_gcta_grm ,plink_ch_gcta_multigrm,plink_ch_gcta_multigrm2,fam_ch_gcta, fam_ch_gcta_mult) { a -> [ a,a,a,a,a[2], a[2]] }
  data_h2gcta1=Channel.fromPath(params.data)
  data_h2gcta1_multi = Channel.fromPath(params.data)
  check_gcta = Channel.create()
  pheno_cols_ch_gcta1=Channel.from(params.pheno.split(",").toList())
  pheno_cols_ch_gcta2=Channel.from(params.pheno.split(",").toList())

  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""

  process  getGctaPhenosCovar {
    input:
      file(covariates) from data_h2gcta1
      file(fam) from fam_ch_gcta
    each this_pheno from pheno_cols_ch_gcta1
    output:
      set val("$this_pheno"),file(phef), file(covfile) into (newdata_ch_gcta_multgrm, newdata_ch_gcta_grm) 
    script:
      base = fam.baseName
      phef = "${this_pheno}_gcta_n.phe"
      covfile = "${this_pheno}_gcta_n.cov"
      """
      all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
                          --pheno ${this_pheno} --phe_out ${phef} --form_out 4 --cov_out $covfile
      """
  }
  if(params.gcta_h2_imp==0){
   //--reml-alg 0
   //Specify the algorithm to run REML iterations, 0 for average information (AI), 1 for Fisher-scoring and 2 for EM. The default option is 0, i.e. AI-REML, if this option is not specified.
   if(params.gcta_h2_grmfile==""){
   process MakeGRM{
     cpus params.gcta_num_cores
     time params.big_time
     memory params.gcta_mem_req
     input:
        set file(bed), file(bim), file(fam)  from plink_ch_gcta_grm
     publishDir "${params.output_dir}/gctagrm", overwrite:true, mode:'copy'
     output :
        set file("tmp.grm.bin"), file("tmp.grm.id"),  file("tmp.grm.N.bin") into gcta_grm
     script :
        plk=bed.baseName
        """
        ${params.gcta_bin} --bfile $plk --make-grm --out tmp --thread-num ${params.gcta_num_cores}
        """
   }
  }else{
  gcta_grm=Channel.fromPath("${params.gcta_h2_grmfile}.grm.id").combine(Channel.fromPath("${params.gcta_h2_grmfile}.grm.bin")).combine(Channel.fromPath("${params.gcta_h2_grmfile}.grm.N.bin"))
  }
   filegrmcta_gcta=gcta_grm.combine(newdata_ch_gcta_grm)
   process doGRLEM_GCTA{
     cpus params.gcta_num_cores
     time params.big_time
     memory params.gcta_mem_req
     input:
        set file(grm1),file(grm2), file(grm3),pheno, file(phef),file(covfile) from filegrmcta_gcta
    publishDir "${params.output_dir}/gcta", overwrite:true, mode:'copy'
     output :
       file("$output"+".hsq")
       file("${output}_gcta.stat") into report_gcta
     script :
        output=pheno.replace('_','-')+"_gcta"
        cov = (params.covariates!="") ?  " --qcovar $covfile " : ""
        """
        ${params.gcta_bin} --reml  --grm tmp --pheno $phef  --grm-cutoff ${params.grm_cutoff}  --thread-num ${params.gcta_num_cores}  --out $output --reml-alg ${params.gcta_reml_alg} $cov
        if [ ! -f $output".hsq" ]
        then
        cat outmultgrlm > $output".hsq"
        fi
        format_correlation.r  ${output}.hsq ${output}_gcta.stat gcta ${pheno} None
        """
  }



  }

   // multi grm in the case
  if(params.gcta_h2_imp==1){
   //case of file with 
   if(params.gcta_h2_grmfile!=""){
    filemultigrmcta_gctai= Channel.fromPath(params.gcta_h2_grmfile)
   }else{
       process GCTAComputeMultiGRM{
          cpus gcta_num_cores
          time params.big_time
          memory params.gcta_mem_reqmgrm
          input :
            set file(bed),file(bim), file(fam) from plink_ch_gcta_multigrm 
          publishDir "${params.output_dir}/gcta/grlem/", overwrite:true, mode:'copy'
          output :
            file("$out"+".score.ld") into grlmscoreld
          script :
            plk=bed.baseName
            out="ldscore_"+params.gcta_h2_ldscore
            """ 
            ${params.gcta_bin} --bfile $plk --ld-score-region ${params.gcta_h2_ldscore} --out $out  --thread-num ${gcta_num_cores}
            """
       } 

       process GCTAStrat {
          input :
             file(scoreld) from grlmscoreld
          publishDir "${params.output_dir}/gcta/grlem/", overwrite:true, mode:'copy'
          output :
             file("snp_group_*") into grlmstroupescore
          script :
            """
            stratify_segment.py --inp_ld  $scoreld --nb_split 4 --out snp_group
            """
       }
       grlmstroupescorec=grlmstroupescore.flatMap().combine(plink_ch_gcta_multigrm2)
       process GCTAGRMByFile{
          input :
            set file(grmfile),file(bed),file(bim),file(fam) from grlmstroupescorec
          publishDir "${params.output_dir}/gcta/grlem/", overwrite:true, mode:'copy'
          output :
             set file("${out}.grm.bin"), file("${out}.grm.id"), file("${out}.grm.N.bin"), file("${out}.log") into resgrlmfile
             val("${out}") into resgrlmhead 
          script :
            plk=bed.baseName
            out=grmfile.baseName
            """
            ${params.gcta_bin} --bfile $plk --extract $grmfile --make-grm --out $out 

            """
       }
      resgrlmfilec=resgrlmfile.collect()
      resgrlmhead=resgrlmhead.collect()
      
      process MergeFile{
          input :
          val(linfo) from resgrlmfilec
          output :
            file("$fileout") into filemultigrmcta_gctai
          script :
            fileout="g"
            lfile2=linfo.flatten().join("\n")
            """
             echo \"\"\"$lfile2 \"\"\"|grep \".log\"|sed 's/.log//g'> $fileout
           """
     } 
   }
   filemultigrmcta_gctai.into{ filemultigrmcta_gcta; filemultigrmcta_gcta_cor}
   filemultigrmcta_gcta=filemultigrmcta_gcta.combine(newdata_ch_gcta_multgrm)
   process doMultiGRM{
     cpus params.gcta_num_cores
     time params.big_time
     memory params.gcta_mem_req
     input:
        set file(listfile),pheno, file(phef),file(covfile) from filemultigrmcta_gcta
    publishDir "${params.output_dir}/gcta", overwrite:true, mode:'copy'
     output :
       file("$output"+".hsq")
     script :
        output=pheno.replace('_','-')+"_gcta"
        covargcta= (params.covariates=="")? "" : " --qcovar  $covfile "
        """
        ${params.gcta_bin} --reml ${params.multigrm_opt} --mgrm $listfile --pheno $phef  --thread-num ${params.gcta_num_cores}  --out $output $covargcta &> outmultgrlm
        if [ ! -f $output".hsq" ]
        then
        cat outmultgrlm > $output".hsq"
        fi
        """
  }
  }
  listpheno=params.pheno.split(",")
  nbpheno=listpheno.size()
  if(params.gcta_h2_multi==1){
    if(nbpheno<2){
        error("\n\n------\npheno number < 2\n\n---\n")
    }  
  process  getGctaPhenosCovarMulti {
    input:
      file(covariates) from data_h2gcta1_multi
      file(fam) from fam_ch_gcta_mult
    output:
      set file(phef), file(covfile) into newdata_ch_gcta_multgrm_multi
    script:
      base = fam.baseName
      phef = "all_gcta_n.phe"
      covfile = "all_gcta_n.cov"
      """
      all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
                          --pheno ${params.pheno} --phe_out ${phef} --form_out 4 --cov_out $covfile
      """
  }
   list_Cor=[]
   for (i = 0; i <(nbpheno-1); i++){
   for( j = i+1; j<nbpheno; j++){
      list_Cor.add([i, j])
   }
   }
   process doMultiGRMCor{
     cpus params.gcta_num_cores
     time params.big_time
     memory params.gcta_mem_req
     input :
       set file(phef), file(covfile) from newdata_ch_gcta_multgrm_multi
       file(filemult) from filemultigrmcta_gcta_cor
    publishDir "${params.output_dir}/gcta", overwrite:true, mode:'copy'
     each poss from list_Cor
     output:
       file("$output"+".hsq")
     script :
        pos=poss[0]
        pos2=poss[1]
        output=listpheno[pos].replace('_','-')+"_"+listpheno[pos2].replace('_','-')
        cov = (params.covariates!="") ?  " --qcovar $covfile " : ""
        pos=pos+1
        pos2=pos2+1
        """
        ${params.gcta_bin} --reml --mgrm $filemult --pheno $phef --thread-num ${params.gcta_num_cores}  --out $output  --reml-bivar $pos $pos2  --reml-alg ${params.gcta_reml_alg}  ${params.gcta_opt_multigrm_cor} $cov  &> outmultgrlm  
        if [ ! -f $output".hsq" ]
        then
        cat outmultgrlm > $output".hsq"
        fi
        """ 
    }

  }

}

//////
///
///
///  GEMMA 
///
//////
  def newNamePheno(Pheno){
      SplP=Pheno.split(',')
      for (i = 0; i <SplP.size(); i++) {
         SplP[i]=(i+1)+"-"+SplP[i]
      }
      return(SplP)
  }

if(params.gemma_h2==1){
    bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
    bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
    fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
    gemma_assoc_ch= Channel.create()
    Channel.from(file(bed),file(bim),file(fam)).buffer(size:3)
        .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
        .set { gemma_assoc_ch }

    rel_ch_gemma= Channel.create()
    gem_ch_gemma = Channel.create()
    gemma_assoc_ch.separate (rel_ch_gemma, gem_ch_gemma) { a -> [ a, a] }

  if(params.gemma_mat_rel==""){
   if(params.file_rs_buildrelat==""){
        filers_matrel_mat_gem=file('NO_FILE')
     }else{
        filers_matrel_mat_gem=Channel.fromPath(params.file_rs_buildrelat)
   }

  process getGemmaRel {
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time params.big_time
    input:
       file plinks from rel_ch_gemma
       file file_rs from filers_matrel_mat_gem
    output:
       file("output/${base}.*XX.txt") into (rel_mat_ch, rel_mat_ch_gxe)
    script:
       base = plinks[0].baseName
       famfile=base+".fam"
       rs_list = params.file_rs_buildrelat!="" ? " -snps $file_rs " : ""
       """
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       cat $famfile |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
       ${params.gemma_bin} -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3 $rs_list
       """
  }
  }else{
   rel_mat_ch=Channel.fromPath(params.gemma_mat_rel, checkIfExists:true)
   rel_mat_ch_gxe=Channel.fromPath(params.gemma_mat_rel, checkIfExists:true)
  }

  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""
  ind_pheno_cols_ch = newNamePheno(params.pheno)
  data_ch = file(params.data)

  gwas_type_gem1=Channel.from(params.gemma_h2_typeest.split(",")).toList()

  process doGemmah2 {
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(covariates) from data_ch
      file(rel) from rel_mat_ch
      set file(bed),file(bim),file(fam) from gem_ch_gemma
    each this_pheno from ind_pheno_cols_ch
    each gemtype from gwas_type_gem1
   publishDir "${params.output_dir}/gemma", overwrite:true, mode:'copy'
    output:
      file("output/${out}.log.txt") 
      file("${out}_gemma.stat") into report_gemma
    script:
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+-/,"")
       ourpheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
       data_nomissing     = "pheno-"+our_pheno+".pheno"
       list_ind_nomissing = "lind-"+our_pheno+".lind"
       rel_matrix         = "newrel-"+our_pheno+".rel"
       gemma_covariate    = "${our_pheno}.gemma_cov"
       phef               = "${our_pheno}_n.phe"
       covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
       out                = "$our_pheno"+"_type"+gemtype
       dir_gemma          =  "gemma"
       base = bed.baseName
       inp_fam = base+".fam"
       newbase= "tmp"
      // check if list phenotype is same order than rel matric
       """
       list_ind_nomissing.py --data $covariates --inp_fam $inp_fam $covariate_option --pheno $ourpheno3 --dataout $data_nomissing \
                             --lindout $list_ind_nomissing
       gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase 
       all_covariate.py --data  $data_nomissing --inp_fam  $newbase".fam" $covariate_option --cov_out $gemma_covariate \
                          --pheno $our_pheno2 --phe_out ${phef} --form_out 1
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       ${params.gemma_bin} ${covar_opt_gemma}  -k $rel_matrix  -n 1 -p $phef -o $out -maf 0.0000001 -vc $gemtype
       format_correlation.r output/$out".log.txt" $out"_gemma.stat" gemma ${our_pheno} $gemtype
       rm $newbase".fam" $newbase".bim" $newbase".bed"
       """

  }
}else{
report_gemma=Channel.empty()
}
if(params.gemma_h2_pval==1){

    bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
    bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
    fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()
    gemmapval_assoc_ch= Channel.create()
    Channel.from(file(bed),file(bim),file(fam)).buffer(size:3)
        .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
        .set { gemmapval_assoc_ch }

gwas_file_gem=Channel.from(params.file_gwas.split(",")).flatMap{it->file(it,checkIfExists:true)}.combine(gemmapval_assoc_ch)

// for 2 we need a zcat file 

typegemmah2=params.gemma_h2_typeest.split(",")
if("2" in typegemmah2){
println "warning for hertitabilies with pvalue option 2 of gemma not implemented issue with wcat option (see manuals)"
}
//gwas_type_gem1=Channel.from(params.gemma_h2_typeest.split(",")).toList()
gwas_type_gem2=Channel.from("1".split(",")).toList()
//gwas_type_gem2=["1"]

process DoGemmah2Pval{
   memory params.gemma_mem_req
   cpus params.gemma_num_cores
   time   params.big_time
   input :
      set file(gwas),file(bed),file(bim),file(fam) from gwas_file_gem
   publishDir "${params.output_dir}/gemmapval", overwrite:true, mode:'copy'
   each gemtype from gwas_type_gem2
   output :
       file("output/*")
       file("${out}_gemmah2.stat") into report_gemmah2
   script :
     NInfo=params.head_n=="" ? " --n_header ${params.head_n}   " : ""
     out=gwas.baseName+"_gemm_"+gemtype.replace('_','-')
     plkbas=bed.baseName
     newplkbas=plkbas+"_new"
     //error! Number of columns in the wcat file does not match that of cat file.error! fail to read files. 
     //WCAT=gemtype=="2" ? " --wcat "
     //This analysis option requires marginal z-scores from the study and individual-level genotypes froma random subset of the study (or a separate reference panel).
     """
     gemma_format_pval.py --inp_asso $gwas --out $gwas".new"  --rs_header ${params.head_rs} --a1_header ${params.head_A1} --a2_header ${params.head_A2} --se_header ${params.head_se} --chro_header ${params.head_chr} --beta_header ${params.head_beta} --bfile $plkbas --threads ${params.gemma_num_cores}
     plink -bfile $plkbas --extract listrs.rs --make-bed  --out $newplkbas --keep-allele-order --threads ${params.gemma_num_cores} 
     cp $newplkbas".fam" $newplkbas".fam.tmp" 
     awk \'{\$6=1;print \$0}\' $newplkbas".fam.tmp" > $newplkbas".fam"
     export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
     ${params.gemma_bin} -beta $gwas".new" -bfile  $newplkbas -vc $gemtype -o $out
     format_correlation.r output/$out".log.txt" $out"_gemmah2.stat" gemmah2 $gwas $gemtype
     rm $newplkbas".fam" $newplkbas".bim" $newplkbas".bed"
     """
}



}

report_ch = report_ldsc.flatten().mix(report_gemma.flatten()).mix(report_bolt.flatten()).mix(report_gcta.flatten()).mix(report_gemmah2.flatten()).toList()

process MergeH2{
   input :
       file(allfile) from report_ch    
   publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
   output :       
       file("{params.output}.*")
   script :
       allfile=allfile.join(',')
       """
       merge_allfile.r $allfile ${params.output}
       """

}

