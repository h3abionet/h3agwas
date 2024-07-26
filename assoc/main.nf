#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Jean-Tristan Brandenburg
 *      Scott Hazelhurst
 *      Shaun Aron
 *   	Rob Clucas
 *      Eugene de Beste
 *      Lerato Magosi
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2020
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

def strmem(val){
 return val as nextflow.util.MemoryUnit
}



//---- General definitions --------------------------------------------------//

import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;
nextflow.enable.dsl = 1 

if (!workflow.resume) {
    def dir = new File(params.output_dir)
    if (dir.exists() && dir.directory && (!(dir.list() as List).empty)) {
       println "\n\n============================================"
       println "Unless you are doing a -resume, the output directory should be empty"
       println "We do not want to overwrite something valuable in "+params.output_dir
       println "Either clean your output directory or check if you meant to do a -resume"
       System.exit(-1)
    }
}


/*definition*/
def errormess(message,exitn=0){
    if(message=="")return(0)
    println(message)
    System.exit(exitn)
}


def checkparams(param, namesparam, type, min=null, max=null, possibleval=null, notpossibleval=null) {
  messageerror=""
  if(param==null){
    messageerror+="error :--"+namesparam+" is null "
  } else {
    if(!(param.getClass() in type)){
   messageerror+="error :--"+namesparam+" must be a "+ type
     if(param.getClass()==Boolean)messageerror+=", but no parameters given"
     else messageerror+=" but type is "+param.getClass()+" value "+ param
   }else{
   if(min && param<min)messageerror+="\nerror : --"+namesparam+" < min value :"+param +" < "+min
   if(max && param>max)messageerror+="\nerror : --"+namesparam +"> maxvalue :" + param+" > "+max
   if(possibleval && !(param in possibleval))messageerror+="\nerro : --"+namesparam +" must be one the value :"+possibleval.join(',')
   }
   }
    errormess(messageerror,2)
}


def checkmultiparam(params, listparams, type, min=null, max=null, possibleval=null, notpossibleval=null){
 messageerror=""
 for(param in listparams){
   if(params.containsKey(param)){
     checkparams(params[param], param, type, min=min, max=max, possibleval=possibleval, notpossibleval=notpossibleval)
   }else{
     messageerror+="param :"+param+" not initialize\n"
   }
 }
 errormess(messageerror, 2)

}








def helps = [ 'help' : 'help' ]

allowed_params_input = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates", "work_dir", "scripts",  "high_ld_regions_fname", "accessKey", "access-key", "secretKey", "secret-key", "region",  "pheno","big_time", "gemma_mat_rel", "file_rs_buildrelat","genetic_map_file", "rs_list",   "gemma_bin", "bgen", "bgen_sample",   "list_bgen", "exclude_snps", "bolt_impute2filelist", "bolt_impute2fidiid", "bolt_otheropt","bolt_bin", "bolt_ld_scores_col" , "bolt_ld_scores_col",  "bolt_impute2filelist", "bolt_impute2fidiid", "bolt_otheropt","bolt_bin", 'gxe','fastlmmc_bin','list_vcf', 'vcf_field', "regenie_otheropt_step1","regenie_otheropt_step2", "gcta_bin", "AMI", "instance-type", "boot-storage-size", "sharedStorageMount", "instanceType", "bolt_ld_score_file" , "saige_impute_method", "gcta_grmfile","fastgwa_type", "snps_include_rel", "sample_snps_rel_paramplkl", "shared-storage-mount", "queue", 'output_testing', "chrom","regenie_bin" , "saige_bin_fitmodel", "saige_bin_spatest", "snp_rel_param_plk", "saige_otheropt", "saige_otheropt_step1", "saige_otheropt_step2"]

allowed_params_input_mp = ["bolt_covariates_type", 'covariates_type']
allowed_params=allowed_params_input
allowed_params+=allowed_params_input_mp
allowed_params_cores=["gemma_num_cores", "max_plink_cores", "bolt_num_cores", 'fastlmm_num_cores', 'saige_num_cores',"regenie_num_cores", "fastgwa_num_cores"]
allowed_params+=allowed_params_cores
allowed_params_intother=["max_forks", "mperm", "regenie_bsize_step1", "regenie_bsize_step2", "grm_nbpart", "thin", "regenie_bsize", "gemma_relopt", "gemma_lmmopt"]
allowed_params+=allowed_params_intother
allowed_params_bolother=["adjust", "mperm", "sample_snps_rel","bolt_use_missing_cov", 'gemma_multi', 'pheno_bin', 'fastlmm_multi', "regenie_loco", "sexinfo_available", "print_pca", "saige_loco","saige_imputed_data"]
allowed_params+=allowed_params_bolother
allowed_params_float=["cut_maf", "bgen_mininfo",  "grm_cutoff", "grm_maf","vcf_minmac", "cut_geno"]
allowed_params+=allowed_params_float
allowed_params_memory=["gemma_mem_req" , "plink_mem_req", "other_mem_req", "bolt_mem_req", 'fastlmm_mem_req', 'saige_mem_req', "regenie_mem_req", "fastgwa_mem_req", "bootStorageSize", "bootStorageSize", "boot-storage-size", "sharedStorageMount", "other_process_mem_req"]
allowed_params+=allowed_params_memory
allowed_params_test=["gemma", "linear","logistic","assoc","fisher",  "saige", "boltlmm", 'gemma_gxe', 'plink_gxe', 'fastlmm', "regenie", "fastgwa", "regenie_gxe","cmh", "model"]
allowed_params+=allowed_params_test








def params_help = new LinkedHashMap(helps)


filescript=file(workflow.scriptFile)
projectdir="${filescript.getParent()}"
dummy_dir="${projectdir}/../qc/input"

params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output_testing = "cleaned"
params.thin       = 0
params.covariates = ""
params.chrom      = ""
params.print_pca = 1
params.file_rs_buildrelat = ""
params.genetic_map_file = ""
params.list_vcf=""
params.vcf_field="DS"
params.vcf_minmac=0.5
outfname = params.output_testing
params.bgen=""
params.bgen_sample=""
params.bgen_mininfo=0.6
params.list_bgen=""
params.snp_rel_param_plk=" --maf 0.01 --mac 10 "


params.regenie_bin="regenie"
params.regenie_bsize=100
params.regenie_bsize_step1=0
params.regenie_bsize_step2=0
params.regenie_otheropt_step1=""
params.regenie_otheropt_step2=""
params.regenie_loco=1
params.regenie_num_cores=6
params.regenie_mem_req="10GB"
//params.regenie_mafstep1=0.01
params.regenie=0
params.saige_imputed_data=1
params.saige_impute_method="best_guess" //best_guess, mean or minor
params.saige_otheropt_step1=""
params.saige_otheropt=""
params.saige_otheropt_step2=""
params.cut_geno = 0.9


/* Defines the path where any scripts to be executed can be found.
 */



/* Do permutation testing -- 0 for none, otherwise give number */
params.mperm = 00

/* Adjust for multiple correcttion */
params.adjust = 0

supported_tests_all = ["assoc","fisher","model","cmh","linear","logistic","boltlmm", "fastlmm", "gemma", "gemma_gxe", 'saige', "regenie"]


params.assoc     = 0
params.fisher   = 0
params.cmh     =  0
params.model   =  0
params.linear   = 0
params.logistic = 0
params.gemma = 0
params.saige=0

params.gemma_multi=0
params.gemma_mem_req = "6GB"
params.gemma_relopt = 1
params.gemma_lmmopt = 4
params.gemma_mat_rel = ""
params.gemma_num_cores = 8
params.pheno = "_notgiven_"

//
params.saige_bin_fitmodel="step1_fitNULLGLMM.R"
params.saige_bin_spatest="step2_SPAtests.R"
params.saige_loco=1
params.saige_mem_req='10GB'
params.saige_num_cores=10

if (params.pheno == "_notgiven_") {
  println "No phenotype given -- set params.pheno";
  System.exit(-2);
}
  

/*JT Append initialisation variable*/
params.bolt_covariates_type = ""
params.bolt_ld_score_file= ""
params.bolt_ld_scores_col=""
params.boltlmm = 0
params.bolt_num_cores=8
params.bolt_mem_req="6GB"
params.bolt_use_missing_cov=0
params.exclude_snps=""
params.bolt_impute2filelist=""
params.bolt_impute2fidiid=""
params.bolt_otheropt=""
/*fastlmm param*/
params.fastlmm = 0
params.fastlmm_num_cores=8
params.fastlmm_mem_req="15GB"
params.fastlmm_multi = 0 

params.fastlmmc_bin ="fastlmmc"
params.bolt_bin ="bolt"
params.gemma_bin ="gemma"

/*gxe param : contains column of gxe*/
params.gemma_gxe=0
params.plink_gxe=0
params.regenie_gxe=0
params.max_plink_cores = 4
params.rs_list=""
params.gxe=""

/**/
params.fastgwa=0
params.fastgwa_mem_req="10G"
params.fastgwa_num_cores=5
params.grm_nbpart=100
params.grm_maf = 0.01
params.gcta_bin = "gcta64"
params.fastgwa_type="--fastGWA-mlm-exact"
params.grm_cutoff =  0.05
params.covariates_type=""
params.gcta_grmfile=""
params.sample_snps_rel=0
params.sample_snps_rel_paramplkl="100 20 0.1"
params.pheno_bin=0
params.snps_include_rel=""


params.input_pat  = 'raw-GWA-data'



params.plink_mem_req = '6GB' // how much plink needs for this
params.other_process_mem_req = '10G' // how much other processed need


plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores 

params.help = false

/*check param*/
params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "\nUnknown parameter : Check parameter <$parm>\n";
  }
}

checkmultiparam(params,allowed_params_input, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_input_mp, [java.lang.String,java.lang.Integer], min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_memory, java.lang.String, min=null, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_cores, java.lang.Integer, min=1, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_intother, java.lang.Integer, min=0, max=null, possibleval=null, notpossibleval=null)
checkmultiparam(params,allowed_params_bolother, java.lang.Integer, min=0, max=null, possibleval=[0,1], notpossibleval=null)
checkmultiparam(params,allowed_params_test, java.lang.Integer, min=0, max=null, possibleval=[0,1], notpossibleval=null)
checkmultiparam(params,allowed_params_float, [java.lang.Float, java.lang.Integer, java.math.BigDecimal, java.lang.Double], min=0, max=null, possibleval=null, notpossibleval=null)

data_ch_pheno = Channel.fromPath(params.data, checkIfExists:true)
data_ch_show = Channel.fromPath(params.data, checkIfExists:true)
data_ch_gemma = Channel.fromPath(params.data, checkIfExists:true)
data_ch_rsselect= Channel.fromPath(params.data, checkIfExists:true)

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


def fileColExists = { fname, pname, cname ->
  if (fname.contains("s3://")){
    println "The file <$fname> is in S3 so we cannot do a pre-check";
    return;
  }
  if (fname.contains("az://")){
    println "The file <$fname> is in Azure so we cannot do a pre-check";
    return;
  }
  f = new File(fname)
  if (! f.exists() && ! fname.contains("s3://") && ! fname.contains("az://")) {
     error("\n\nThe file <${fname}> given for <${pname}> does not exist")
    } else {
      def line  
      f.withReader { line = it.readLine() }  
      // now get the column headers
      fields = line.split()
      // now separate the column
      cols = cname.split(",")
      cols.each { col -> 
	det = col.split("/")
	if ((det[0].length()>0) && (! fields.contains(det[0])))
	  error("\n\nThe file <${fname}> given for <$pname> does not have a column <${det}>\n")
      }
    }
}

fileColExists(params.data,"${params.data} - covariates", params.covariates)
fileColExists(params.data,"${params.data} - phenotypes", params.pheno)
fileColExists(params.data,"${params.data} - gxe", params.gxe)

covs =  params.covariates.split(",")
params.pheno.split(",").each { p ->
  if (covs.contains(p)) {
    println("\n\nThe phenotype <$p> is also given as a covariate -- this seems like a very bad idea")
    sleep(10)
  }
}




//---- Modification of variables for pipeline -------------------------------//


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




bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString().replaceFirst(/^az:/, "az:/").replaceFirst(/^s3:/, "s3:/")


gemma_assoc_ch = Channel.create()
/*JT initatilisation of boltlmm_assoc_ch*/
boltlmm_assoc_ch = Channel.create()
fastlmm_assoc_ch = Channel.create()
rel_ch_fastlmm = Channel.create()

pca_in_ch = Channel.create()
assoc_ch  = Channel.create()
assoc_ch_gxe  = Channel.create()
assoc_ch_gxe_freq=Channel.create()
grlm_assoc_ch = Channel.create()
ch_bolt_snpchoice=Channel.create()
fastgwa_assoc_ch = Channel.create()
raw_src_ch= Channel.create()
ch_select_rs_format=Channel.create()
ch_format_ldscore=Channel.create()
ch_saige_heritability=Channel.create()
ch_saige_assoc=Channel.create()
ch_regenie_assoc=Channel.create()

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


println "\nTesting data            : ${params.input_pat}\n"
println "Testing for phenotypes  : ${params.pheno}\n"
println "Using covariates        : ${params.covariates}\n\n"

for(test in allowed_params_test){
 if (params[test]) println "Doing $test testing"
}

if (params.thin)
   thin = "--thin ${params.thin}"
else 
   thin = ""

if (params.chrom) 
   chrom = "--chr ${params.chrom}"
else
   chrom = ""

if (thin+chrom) {
  process thin {
    input: 
      set file(bed), file(bim), file(fam) from raw_src_ch
    output:
      /*JT Append initialisation boltlmm_assoc_ch */
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  ( pca_in_ch, assoc_ch, gemma_assoc_ch, boltlmm_assoc_ch,fastlmm_assoc_ch, rel_ch_fastlmm,assoc_ch_gxe_freq, assoc_ch_gxe, grlm_assoc_ch, fastgwa_assoc_ch, ch_bolt_snpchoice, ch_select_rs_format, ch_saige_heritability, ch_saige_assoc,ch_regenie_assoc)
    script:
       base = bed.baseName
       out  = base+"_t"
       "plink --keep-allele-order --bfile $base $thin $chrom --make-bed --out $out --allow-extra-chr --geno ${params.cut_geno} "
  }

  println "\nData has been thinned or only some chromosomes used  (is the run for test purposes only?)\n"
   


} else {
    /*JT : append boltlmm_assoc_ch and a]*/
    raw_src_ch.separate( pca_in_ch, assoc_ch, assoc_ch_gxe,assoc_ch_gxe_freq,gemma_assoc_ch, boltlmm_assoc_ch, fastlmm_assoc_ch,rel_ch_fastlmm, grlm_assoc_ch,fastgwa_assoc_ch,ch_bolt_snpchoice,ch_select_rs_format, ch_format_ldscore,ch_saige_heritability, ch_saige_assoc,ch_regenie_assoc) { a -> [a,a,a,a,a,a,a,a,a,a,a, a,a,a,a, a] }
}




if(params.print_pca!=0){
   process computePCA {
     cpus max_plink_cores
     memory plink_mem_req
     time   params.big_time
     memory { strmem(plink_mem_req) + 5.GB * (task.attempt -1) }
     errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
     input:
       set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from pca_in_ch
     publishDir params.output_dir, overwrite:true, mode:'copy'
     output:
       set file("${outfname}.eigenval"), file("${outfname}.eigenvec")  \
	    into pca_out_ch
     script:
	 base = "cleaned"
	 prune= "${base}-prune"
	"""
	plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check --maf 0.01  --autosome --allow-extra-chr
	plink --keep-allele-order --bfile ${base} --extract check.prune.in --make-bed --out $prune  --autosome --allow-extra-chr
	plink --threads $max_plink_cores --bfile $prune --pca --out ${outfname}   --autosome --allow-extra-chr
	"""
   }

process drawPCA {
    input:
      set file(eigvals), file(eigvecs) from pca_out_ch
    output:
      set file (output), file ("B040-pca.tex") into report_pca_ch
    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.pdf"
    script:
      base=eigvals.baseName
      cc_fname = 0
      cc       = 0
      col      = 0
      // also relies on "col" defined above
      output="${base}-pca.pdf"
      template "drawPCA.py"
   }
}
else {
  report_pca_ch = Channel.empty()
}




//scott, why are you limited at one core when you don't have permutation?
//num_assoc_cores = params.mperm == 0 ? 1 : Math.min(10,params.max_plink_cores)
num_assoc_cores = params.max_plink_cores
supported_tests = ["assoc","fisher","model","cmh","linear","logistic"]

requested_tests = supported_tests.findAll { entry -> params.get(entry) }


covariate = ""
gotcovar  = 0
pheno     = ""

/*Case where we sample builrelatdness*/
balise_filers_rel=1
if(params.boltlmm+params.gemma+params.fastlmm+params.fastgwa+params.saige+params.gemma_gxe+params.saige+params.regenie>0){
 if(params.file_rs_buildrelat=="" && params.sample_snps_rel==1){
  if(params.snps_include_rel!='') snpinclrel=channel.fromPath(params.snps_include_rel)
  else snpinclrel=channel.fromPath("${dummy_dir}/00")
  process select_rs_format{
     cpus max_plink_cores
     memory plink_mem_req
     time   params.big_time
     input :
       tuple path(bed),path(bim), path(fam) from ch_select_rs_format
       path(data) from  data_ch_rsselect
       path(snp_inclrel) from  snpinclrel
    output:
       file("${prune}.prune.in") into  filers_matrel_mat_fast, filers_matrel_mat_GWA, filers_matrel_mat_gem, filers_matrel_bolt, filers_count_line, filers_her_saige,filers_matrel_regenie
     script:
        base = bed.baseName
        prune= "${base}-prune"
        extract=(params.snps_include_rel=='')? "" : " -extract range $snp_inclrel "
        """
        plink --bfile ${base} --indep-pairwise ${params.sample_snps_rel_paramplkl} --out $prune   --threads ${params.max_plink_cores} $extract  --keep $data  ${params.snp_rel_param_plk} --allow-extra-chr
        """
   }
   //BoltNbMaxSnps=filers_count_line.countLines()
 }else{
/* n*/
 if(params.file_rs_buildrelat==""){
   balise_filers_rel=0
   filers_matrel_mat_fast=Channel.fromPath("${dummy_dir}/0",checkIfExists:true) 
   filers_matrel_mat_GWA=Channel.fromPath("${dummy_dir}/0") 
   filers_matrel_mat_gem=channel.fromPath("${dummy_dir}/0") 
   filers_her_saige=channel.fromPath("${dummy_dir}/00")
   filers_matrel_regenie=channel.fromPath("${dummy_dir}/00")
   if(params.boltlmm==1){
      //BoltNbMaxSnps=1000000
      process buildBoltFileSnpRel{
         memory params.bolt_mem_req
         time   params.big_time
         input:
           set file(bed),file(plinksbim), file(fam) from ch_bolt_snpchoice
         output :
           file(output) into filers_matrel_bolt
         script :
           output=plinksbim.baseName+".rs.choice"
           """
           shuf -n 950000 $plinksbim | awk '{print \$2}' > $output
           """
      }

  }

  }else{
        filers_matrel_mat_fast=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
        filers_matrel_bolt=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
        //if(params.boltlmm==1)BoltNbMaxSnps=CountLinesFile(params.file_rs_buildrelat)
        filers_matrel_mat_GWA=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
        filers_matrel_mat_gem=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
        filers_her_saige=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
        filers_matrel_regenie=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
  }
 }
}

/*format and prepared bgen sample*/
if(params.bgen!="" || params.list_bgen!=''){
    if(params.bgen_sample==''){
      println "params.bgen_sample not initialise when bgen params initial";
      System.exit(-2);
     }
    data_ch_bgen=channel.fromPath(params.data, checkIfExists:true)
    bgensample_ch_i = Channel.fromPath(params.bgen_sample, checkIfExists:true)
    process bgen_formatsample {
        label 'R'
        input :
         path(data) from data_ch_bgen
         path(bgen_sample) from  bgensample_ch_i
       output :
          path(bgen_sample2)  into (bgensample_ch,bgensample_ch2, bgensample_regenie_ch, bgensample_ch_regenie)
          path(bgen_samplesaige2) into  (bgensample_ch_fastgwa,bgensample_ch_saige)
       script :
           bgen_sample2=bgen_sample+".modif"
           bgen_samplesaige2=bgen_sample+"_saige.modif"
           """
           format_samplebgen.r --sample $bgen_sample  --out_sample $bgen_sample2 --out_samplesaige $bgen_samplesaige2 --data $data
           """
    } 

}else{
  bgensample_ch= Channel.fromPath("${dummy_dir}/07", checkIfExists:true)
  bgensample_ch2= Channel.fromPath("${dummy_dir}/07", checkIfExists:true)
  bgensample_ch_fastgwa= Channel.fromPath("${dummy_dir}/07", checkIfExists:true)
  bgensample_ch_saige= Channel.fromPath("${dummy_dir}/07", checkIfExists:true)
  bgensample_ch_regenie= Channel.fromPath("${dummy_dir}/07", checkIfExists:true)


}
 
if (params.data != "") {

   checker(file(params.data))

   if (params.covariates != "") {
      gotcovar = 1
  }


  
  process extractPheno {
    input:
     file(data) from data_ch_pheno
    output:
     file(phenof) into pheno_ch
    script:
     phenof = "pheno.phe"
     all_phenos = params.covariates.length()>0 ? params.pheno+","+params.covariates : params.pheno
     """
     extractPheno.py $data ${all_phenos} $phenof
     """
  }


  pheno_label_ch = Channel.from(params.pheno.split(","))

  process showPhenoDistrib {
    // not sure difference between container and label
    input:
    file(data) from data_ch_show
    output:
      file ("B050*") into pheno_report_ch
    script:
      "phe_distrib.py --pheno ${params.pheno} $data B050 "
  }
}  else {
  pheno_report_ch = Channel.empty()
  pheno_label = ""
  pheno_label_ch = Channel.from("")
}

/*JT : Case fastlmm => if yes*/
if (params.fastlmm == 1) {
  data_ch_fastlmm = Channel.fromPath(params.data, checkIfExists:true)
  if(params.fastlmmc_bin=="")fastlmmc="fastlmmc"
  else fastlmmc=params.fastlmmc_bin

  fam_ch_fast = Channel.create()
  gem_ch_fast2 = Channel.create()
  gem_ch_fast =Channel.create()
  bim_ch_fast_fas = Channel.create()
  fastlmm_assoc_ch.separate (gem_ch_fast2,gem_ch_fast,bim_ch_fast_fas,fam_ch_fast) { a -> [a,a, a[1],a[2]] }


  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""

  process  getFastLmmPhenosCovar {
    input:
      file(covariates) from data_ch_fastlmm
      file(fam) from fam_ch_fast
    output:
      set file(phef), file(covfile) into fastlmm_data_ch
      stdout into pheno_cols_ch_fastlmm
    script:
      base = fam.baseName
      phef = "${base}_fastlmm_n.phe"
      covfile = "${base}_fastlmm_n.cov"
      """
      all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
                          --pheno ${params.pheno} --phe_out ${phef}  --cov_out $covfile --form_out 3
      """
  }

  ind_pheno_cols_ch = Channel.create()
  check = Channel.create()
  pheno_cols_ch_fastlmm.flatMap { list_str -> list_str.split() }.tap ( check) .set { ind_pheno_cols_ch }

  if(params.fastlmm_multi==1){

     process getRelForFastLMM {
        label 'gemma'
	cpus params.fastlmm_num_cores
	memory params.fastlmm_mem_req
	time params.big_time
	input:
	   file plinks from rel_ch_fastlmm
           file file_rs from filers_matrel_mat_fast
	output:
	   file("output/${base}.*XX.txt")
	   file("${rel_fastlmm}") into rel_mat_ch_fastlmm
	script:
	   base = plinks[0].baseName
	   fam = plinks[2]
	   rel_fastlmm="rel_fastlmm.txt"
           rs_list = balise_filers_rel== 1 ? " -snps $file_rs " : ""
	   """
           cat $fam |awk '{print \$1"\t"\$2"\t"0.2}' > pheno
	   export OPENBLAS_NUM_THREADS=${params.fastlmm_num_cores}
	   ${params.gemma_bin} -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3 $rs_list
	   cvt_rel_gemma_fastlmm.py $fam output/${base}.*XX.txt $rel_fastlmm
	   """
	 }


     process getListeChro{
	input :
	  file(BimFile) from bim_ch_fast_fas
	output :
	  stdout into (chrolist,chrolist2)
	script:
	 """
	 cat $BimFile|awk '{print \$1}'|uniq|sort|uniq
	"""
     }

     check2 = Channel.create()
     ListeChro2=chrolist.flatMap { list_str -> list_str.split() }.tap ( check2)


     process doFastlmmMulti{
       label 'fastlmm'
       cpus params.fastlmm_num_cores
       memory params.fastlmm_mem_req
       time   params.big_time
       maxForks params.max_forks
       input:
	 set file (phef), file(covariate) from fastlmm_data_ch
	 file(rel) from rel_mat_ch_fastlmm
	 file(plinks) from  gem_ch_fast
       each this_pheno from ind_pheno_cols_ch
       each chro from ListeChro2
       output:
	 set (our_pheno, file("$out"), val(base)) into (fastlmm_manhatten_chro,fastlmm_manhatten_chro2)
       script:
	 base = plinks[0].baseName
         our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
         our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")

	 covar_opt_fast =  (params.covariates) ?  " -covar newcov.out" : ""
	 newbase=base+"-"+chro
	 out = "$base-$our_pheno"+"-"+chro+".stat"
	 """
	 this_pheno_col=`echo ${this_pheno} | awk -F"@" '{print \$1}'`
         fastlmm_relselind.py --rel $rel --phenofile $phef --relout rel_fastlmm_filter.txt --phenofileout newpheno.out --pospheno \$this_pheno_col --covfile $covariate --covfileout newcov.out
         plink --keep-allele-order --bfile $base --chr $chro --make-bed --out $newbase --threads ${params.fastlmm_num_cores} --keep newpheno.out
	 $fastlmmc -REML -simType RRM -verboseOut -sim $rel -bfile $newbase -pheno ${phef} -simLearnType Full -out $out -maxThreads ${params.fastlmm_num_cores} \
	          $covar_opt_fast  
	 """
       }

     fastlmm_manhatten_chroM=fastlmm_manhatten_chro.groupTuple()
     fastlmm_manhatten_chroM1=fastlmm_manhatten_chro2.groupTuple()


     process doMergeFastlmm{
          input :
	    set (val(this_pheno),list_file, base_list) from fastlmm_manhatten_chroM
	    /* with uniq channels vs 2 => problems*/
	    /*file(plinks) from  gem_ch_fast2*/
	 publishDir "${params.output_dir}/fastlmm", overwrite:true, mode:'copy'
	 output :
	     set val(base), val(our_pheno2), file("$out") into fastlmm_manhatten_ch
	 script :
	     base=base_list[0]
	     //our_pheno = this_pheno.replace(/_|\/np.\w+/,"-").replace(/-$/,"")
	     //our_pheno2 = this_pheno.replace(/_|\/np.\w+/,"-").replace(/-$/,"").replaceAll(/^[0-9]+-/,"")
             our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
             our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")

	     out = "$base-${our_pheno}.stat"
	     fnames = list_file.join(" ")
	     file1  = list_file[0]
	     """
	     head -1 $file1 > $out
	     cat $fnames | grep -v "Chromosome" >> $out
	     """
     }
  }  else { // if not   doing fastlmm_multi

     process doFastlmm{
       maxForks params.max_forks
       label 'fastlmm'
       cpus params.fastlmm_num_cores
       time   params.big_time
       memory params.fastlmm_mem_req
       input:
	 set file(phef), file (covariate) from fastlmm_data_ch
	 file(plinks) from  gem_ch_fast
       publishDir "${params.output_dir}/fastlmm", overwrite:true, mode:'copy'
       each this_pheno from ind_pheno_cols_ch
       output:
         file(out)
	 set val(base), val(our_pheno2), file("$out") into fastlmm_manhatten_ch
       script:
	 base = plinks[0].baseName
//	 our_pheno = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
//	 our_pheno2 = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"").replaceAll(/^[0-9]+-/,"")

         our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
         our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
	 covar_opt_fast =  (params.covariates) ?  " -covar $covariate" : ""
	 out = "$base-$our_pheno"+".stat"
	 """
	 this_pheno_col=`echo ${this_pheno} | awk -F"@" '{print \$1}'`
	 $fastlmmc -REML -simType RRM -verboseOut -bfile $base -pheno ${phef} -simLearnType Full -out $out -maxThreads $params.fastlmm_num_cores \
	           $covar_opt_fast  -mpheno \${this_pheno_col} -bfileSim $base
	 """
     }
  }

  // this part is plotting done for any fastlmm mode
  //, overwrite:true, mode:'copy'
  process showFastLmmManhattan {
    memory params.other_process_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      set val(base), val(this_pheno), file(assoc) from fastlmm_manhatten_ch
    output:
      file("${out}*")  into report_fastlmm_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C051-fastlmm-"+this_pheno
      """
      general_man.py  --inp $assoc --phenoname $this_pheno --out ${out} --chro_header Chromosome --pos_header Position --rs_header SNP --pval_header Pvalue --beta_header SNPWeight --info_prog FastLmm
      """
  }




  // End of FASTLMM
}  else {
  report_fastlmm_ch=Channel.empty()
} 


/*JT : Case boltlmm => if yes*/


   /*JT Fonction to transforme argument for cofactor in gemma
   @Input 
   args: cofactor args separate by a comma
   infoargs: type of cofactor separate by a comma : 0 for qualitative, 1 for quantitative
   output : cofactor for boltlmm was formating to take account qualitative and quantitative
   */
   def boltlmmCofact(args,infoargs) {
      //Method code
      splargs=args.split(",")
      infoargs=""+infoargs
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




if (params.boltlmm == 1) {

  plink_ch_bolt = Channel.create()
  fam_ch_bolt = Channel.create()
  bim_ch_bolt = Channel.create()
  boltlmm_assoc_ch.separate (plink_ch_bolt, fam_ch_bolt, bim_ch_bolt) { a -> [ a, a[2], a[1]] }
  data_ch_bolt = Channel.fromPath(params.data, checkIfExists:true)
  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""

  process  getBoltPhenosCovar {
    input:
      file(covariates) from data_ch_bolt
      file(fam) from fam_ch_bolt
    output:
      file(phef) into newdata_ch_bolt
      stdout into pheno_cols_ch_bolt
    script:
      base = fam.baseName
      phef = "${base}_fastlmm_n.phe"
      """
      all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
                          --pheno ${params.pheno} --phe_out ${phef} --form_out 2
      """
  }

  ind_pheno_cols_ch_bolt = Channel.create()
  check_bolt = Channel.create()
  pheno_cols_ch_bolt.flatMap { list_str -> list_str.split() }.tap ( check_bolt) .set { ind_pheno_cols_ch_bolt }


  if(params.bolt_covariates_type=="" & params.covariates_type!=""){
    bolt_covariates_type=params.covariates_type
  }else{
    bolt_covariates_type=params.bolt_covariates_type
  }
   if (params.covariates) 
      cov_bolt = boltlmmCofact(params.covariates,bolt_covariates_type)
   else
      cov_bolt= ""

   missing_cov=""
   if(params.bolt_use_missing_cov==1)
     missing_cov=" --covarUseMissingIndic "

  pval_head = "P_BOLT_LMM"
  type_lmm="--lmm"

  if(params.exclude_snps)rs_ch_exclude_bolt=Channel.fromPath(params.exclude_snps, checkIfExists:true)
  else rs_ch_exclude_bolt=Channel.fromPath("${dummy_dir}/01", checkIfExists:true) 

  if(params.bolt_impute2filelist!=""){
  Impute2FileList=Channel.fromPath(params.bolt_impute2filelist, checkIfExists:true)
  Impute2FID = Channel.fromPath(params.bolt_impute2fidiid, checkIfExists:true)
  }else{
  Impute2FileList=Channel.fromPath("${dummy_dir}/02", checkIfExists:true)
  Impute2FID = Channel.fromPath("${dummy_dir}/03", checkIfExists:true)
  }
  if(params.bolt_ld_score_file!=""){
     Bolt_ld_scoreI= Channel.fromPath(params.bolt_ld_score_file, checkIfExists:true)
     process format_genetic_ldscore{
       input :
        file(ldscore) from Bolt_ld_scoreI
        set file(bed), file(bim), file(fam) from ch_format_ldscore
       output :
         file(newldscore2) into Bolt_ld_score 
       script :
         plkhead=bed.baseName
         newldscore2=ldscore.baseName+"_format.gz"
         """
         zcat $ldscore| awk '{print \$2"\\t"\$3"\\t"\$3"\\t"\$1}' > $ldscore".tmp.bed" 
         plink -bfile $plkhead --extract range $ldscore".tmp.bed" -out $plkhead".tmp" --make-bed
         ldscore_format.py --ldscore $ldscore --plk_bim $plkhead".tmp.bim" --out $newldscore2
         """

     }
  }else{
     Bolt_ld_score = Channel.fromPath("${dummy_dir}/04",checkIfExists:true) 
  }
  if(params.genetic_map_file!=""){
     Bolt_genetic_map= Channel.fromPath(params.genetic_map_file, checkIfExists:true)
  }else{
     Bolt_genetic_map = Channel.fromPath("${dummy_dir}/05",checkIfExists:true) 
  }
  if(params.bgen!="")bgen_ch=Channel.fromPath(params.bgen, checkIfExists:true)
  else bgen_ch=Channel.fromPath("${dummy_dir}/08", checkIfExists:true)

  if(params.list_bgen!="")bgen_ch=Channel.fromPath(file(params.list_bgen).readLines(), checkIfExists:true).collect()
  process doBoltmm{
    maxForks params.max_forks
    label 'bolt'
    cpus params.bolt_num_cores
    memory params.bolt_mem_req
    time   params.big_time
    input:
      tuple path(plinksbed), path(plinksbim), path(plinksfam) from plink_ch_bolt
      path(phef) from newdata_ch_bolt
      file(rs_exclude) from rs_ch_exclude_bolt
      file(SnpChoiceMod) from filers_matrel_bolt
      file(imp2_filelist) from Impute2FileList 
      file(imp2_fid) from Impute2FID
      file(bolt_ld_score) from Bolt_ld_score
      file(bolt_genetic_map) from Bolt_genetic_map
      path(bgen) from bgen_ch
      path(bgensample) from bgensample_ch
    publishDir "${params.output_dir}/boltlmm", overwrite:true, mode:'copy'
    each this_pheno from ind_pheno_cols_ch_bolt
    output:
      file(outbolt)
      set val(base), val(our_pheno), file("$outf") into bolt_manhatten_ch
      path("$base-${our_pheno2}.*")
    script:
      base = plinksbed.baseName
      our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
      our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
      our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
      outimp1  = (params.bolt_impute2filelist!="") ? "$base-${our_pheno2}.imp.stat" : "$base-${our_pheno2}.stat"
      outimp2  = (params.bgen!="") ? "$base-${our_pheno2}.bgen.stat" : "$base-${our_pheno2}.stat"
      outimp2  = (params.list_bgen!="") ? "$base-${our_pheno2}.bgen.stat" : outimp2
      outbolt     = "$base-${our_pheno2}.stat" 
      outf    = (params.bolt_impute2filelist!="") ? outimp1 : outbolt
      outf    = (params.bgen!="") ? outimp2 : outf
      outf    = (params.list_bgen!="") ? outimp2 : outf
      outReml = "$base-$our_pheno2"+".reml"
      covar_file_bolt =  (params.covariates) ?  " --covarFile ${phef} " : ""
      model_snp  = "--modelSnps=$SnpChoiceMod "
      ld_score_cmd = (params.bolt_ld_score_file!="") ? "--LDscoresFile=$bolt_ld_score" :" --LDscoresUseChip "
      ld_score_cmd = (params.bolt_ld_score_file!="" & params.bolt_ld_scores_col!="") ? "$ld_score_cmd --LDscoresCol=${params.bolt_ld_scores_col}" :" $ld_score_cmd "
      exclude_snp = (params.exclude_snps!="") ? " --exclude $rs_exclude " : ""
      boltimpute = (params.bolt_impute2filelist!="") ? " --impute2FileList $imp2_filelist --impute2FidIidFile $imp2_fid --statsFileImpute2Snps $outimp1 --impute2MinMAF ${params.cut_maf} " : ""
      bgen= (params.bgen!="") ? " --bgenFile $bgen --sampleFile $bgensample  --bgenMinINFO ${params.bgen_mininfo}  --bgenMinMAF ${params.cut_maf} --statsFileBgenSnps $outimp2 " : ""
      listbgen="list_bgen"
      bgen= (params.list_bgen!="") ? " --bgenSampleFileList $listbgen --bgenMinINFO ${params.bgen_mininfo}  --bgenMinMAF ${params.cut_maf} --statsFileBgenSnps $outimp2 " : " $bgen "
      baliselistbgen=(params.list_bgen!="")  ? 1 : 0
      geneticmap = (params.genetic_map_file!="") ?  " --geneticMapFile=$bolt_genetic_map " : ""
      """
      if [ "$baliselistbgen" -eq "1" ]
      then
      ls *.bgen | awk -v sample=${bgensample} '{print \$1" "sample}' > $listbgen
      fi
      BoltNbMaxSnps=`cat  ${SnpChoiceMod}|wc -l`
      bolt.py ${params.bolt_bin} $type_lmm --bfile=$base  --phenoFile=${phef} --phenoCol=${our_pheno3}  --numThreads=$params.bolt_num_cores $cov_bolt $covar_file_bolt --statsFile=$outbolt $missing_cov $bgen $boltimpute $ld_score_cmd $exclude_snp  $model_snp $geneticmap  ${params.bolt_otheropt} --maxModelSnps=\$BoltNbMaxSnps --lmmForceNonInf
      cp tmp.stderr $base-${our_pheno2}.boltstderr
      cp .command.out $base-${our_pheno2}.boltstdout
      """
  }

  process showBoltmmManhattan {
   memory params.other_process_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      set val(base), val(this_pheno), file(assoc) from bolt_manhatten_ch
    output:
      file("${out}*")  into report_bolt_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C052-boltlmmm-"+this_pheno
      """
      general_man.py  --inp $assoc --phenoname $this_pheno --out ${out} --chro_header CHR --pos_header BP --rs_header SNP --pval_header $pval_head --beta_header BETA --info_prog BoltLMM
      """
  }


}else {
  report_bolt_ch = Channel.empty()
}
 

  def newNamePheno(Pheno){
      SplP=Pheno.split(',')
      for (i = 0; i <SplP.size(); i++) {
         SplP[i]=(i+1)+"@@@"+SplP[i]
      }
      return(SplP)
  }



if (params.gemma+params.gemma_gxe>0) {

  rel_ch_gemma = Channel.create()
  gem_ch_gemma = Channel.create()
  bim_ch_fast_gem = Channel.create()
  gem_ch_gemma_gxe = Channel.create()
  gemma_assoc_ch.separate (rel_ch_gemma, gem_ch_gemma, gem_ch_gemma_gxe, bim_ch_fast_gem) { a -> [a, a, a,a[1]] }
  if(params.gemma_mat_rel==""){
     process getGemmaRel {
       label 'gemma'
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time params.big_time
       input:
	  file plinks from rel_ch_gemma
	  file file_rs from filers_matrel_mat_gem
       publishDir "${params.output_dir}/gemma/rel", overwrite:true, mode:'copy'
       output:
	  file("output/${base}.*XX.txt") into (rel_mat_ch, rel_mat_ch_gxe)
       script:
	  base = plinks[0].baseName
	  famfile=base+".fam"
	  rs_list = balise_filers_rel==1 ? " -snps $file_rs " : ""
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

   }

   if(params.gemma+params.gemma_gxe>0 & params.gemma_multi==1){
     process getListeChroGem{
        input :
          file(BimFile) from bim_ch_fast_gem
        output :
          stdout into (chrolist_gem,chrolist2_gem, chrolisti_gem_gxe)
        script:
         """
         cat $BimFile|awk '{print \$1}'|uniq|sort|uniq
        """
      }


      check2 = Channel.create()
      list_chro_gemma=chrolist_gem.flatMap { list_str -> list_str.split() }.tap ( check2)
      check2 = Channel.create()
      list_chro_gemma_gxe=chrolisti_gem_gxe.flatMap { list_str -> list_str.split() }.tap ( check2)



   }

   
  if (params.gemma==1 & params.gemma == 1){
      ind_pheno_cols_ch = newNamePheno(params.pheno)

      if (params.covariates)
        covariate_option = "--cov_list ${params.covariates}"
      else
         covariate_option = ""
      if(params.rs_list=="")
        rsfile=Channel.fromPath("${dummy_dir}/06", checkIfExists:true)
       else
        rsfile=Channel.fromPath(params.rs_list, checkIfExists:true)
 


      if(params.gemma_multi==1){

      check2 = Channel.create()

  process doGemmaChro{
       maxForks params.max_forks
       label 'gemma'
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time   params.big_time
       input:
	 file(covariates) from data_ch_gemma
	 file(rel) from rel_mat_ch
	 file(plinks) from  gem_ch_gemma
	 file(rsfilelist) from rsfile
       each this_pheno from ind_pheno_cols_ch
       each chro from list_chro_gemma
       output:
	 file("${dir_gemma}/${out}.log.txt")
	 set  val(our_pheno),file("${dir_gemma}/${out}.assoc.txt"), val(base) into gemma_manhatten_ch_chro
       script:
	  our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
	  our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
	  our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
	  data_nomissing     = "pheno-"+our_pheno+".pheno"
	  list_ind_nomissing = "lind-"+our_pheno+".lind"
	  rel_matrix         = "newrel-"+our_pheno+".rel"
	  base               =  plinks[0].baseName
	  inp_fam            =  base+".fam"
	  newbase            =  base+"-"+our_pheno
	  newfam             =  newbase+".fam"
	  gemma_covariate    = "${newbase}.gemma_cov"
	  phef               = "${newbase}_n.phe"
	  covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
	  rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
	  out                = "$base-$our_pheno-$chro"
	  dir_gemma          =  "gemma"
	  """
	  hostname
	  list_ind_nomissing.py --data $covariates --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout $data_nomissing \
				--lindout $list_ind_nomissing
	  gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
	  plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem} --chr $chro
	  all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
			     --pheno $our_pheno2 --phe_out ${phef} --form_out 1
	  export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
	  ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf ${params.cut_maf}
	  mv output ${dir_gemma}
	  rm $rel_matrix
	  rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
	  """
     }


     gemma_manhatten_ch_chro_merge=gemma_manhatten_ch_chro.groupTuple()

     process doMergeGemma{
	     input :
	       set (val(this_pheno),file(list_file), base_list) from  gemma_manhatten_ch_chro_merge
	    publishDir "${params.output_dir}/gemma", overwrite:true, mode:'copy'
	    output :
		set val(base), val(our_pheno2), file("$out") into gemma_manhatten_ch
	    script :
		base=base_list[0]
		our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
		our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
		out = "$base-${our_pheno}.gemma"
		fnames = list_file.join(" ")
		file1  = list_file[0]
		"""
		head -1 $file1 > $out
		cat $fnames | grep -v "p_wald" >> $out
		"""
	}

     } else {

     process doGemma{
       maxForks params.max_forks
       label 'gemma'
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time   params.big_time
       input:
	 file(covariates) from data_ch_gemma
	 file(rel) from rel_mat_ch
	 file(plinks) from  gem_ch_gemma
	 file(rsfilelist) from rsfile
       each this_pheno from ind_pheno_cols_ch
       publishDir params.output_dir, overwrite:true, mode:'copy'
       output:
	 file("${dir_gemma}/${out}.log.txt")
	 set val(newbase), val(our_pheno), file("${dir_gemma}/${out}.assoc.txt") into gemma_manhatten_ch
       script:
	  our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
	  our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
	  our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
	  data_nomissing     = "pheno-"+our_pheno+".pheno"
	  list_ind_nomissing = "lind-"+our_pheno+".lind"
	  rel_matrix         = "newrel-"+our_pheno+".rel"
	  base               =  plinks[0].baseName
	  inp_fam            =  base+".fam"
	  newbase            =  base+"-"+our_pheno
	  newfam             =  newbase+".fam"
	  gemma_covariate    = "${newbase}.gemma_cov"
	  phef               = "${newbase}_n.phe"
	  covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
	  rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
	  out                = "$base-$our_pheno"
	  dir_gemma          =  "gemma"
	  """
	  list_ind_nomissing.py --data $covariates --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout $data_nomissing \
				--lindout $list_ind_nomissing
	  gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
	  plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem}
	  all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
			     --pheno $our_pheno2 --phe_out ${phef} --form_out 1
	  export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
	  ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf ${params.cut_maf}
	  mv output ${dir_gemma}
	  rm $rel_matrix
	  rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
	  """
      }
   }


   process showGemmaManhattan {
    memory params.other_process_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    label 'bigMem'
    input:
      set val(base), val(this_pheno), file(assoc) from gemma_manhatten_ch
    output:
      file("${out}*")  into report_gemma_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C053$this_pheno"
      """
      gemma_man.py  $assoc $this_pheno ${out}
      """
  }


  } else {
     report_gemma_ch = Channel.empty()
  }


  if (params.gemma_gxe == 1){
     data_ch_gxe = Channel.fromPath(params.data, checkIfExists:true)
   
  if (params.gemma_gxe) 
    gxe_option = "--gxe ${params.gxe}"
  else 
    gxe_option = ""
   if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""
   if(params.rs_list==""){
        rsfile=Channel.fromPath("${dummy_dir}/07",checkIfExists:true) 
     }else{
        rsfile=Channel.fromPath(params.rs_list,checkIfExists:true)
   }

   if(params.gemma_multi==1){
     ind_pheno_cols_ch_gxe_multi = newNamePheno(params.pheno)

    process doGemmaGxEChro{
       maxForks params.max_forks
       cpus params.gemma_num_cores
       memory params.gemma_mem_req
       time   params.big_time
       input:
	 file(covariates) from data_ch_gxe
	 file(rel) from rel_mat_ch_gxe
	 file(plinks) from  gem_ch_gemma_gxe
	 file(rsfilelist) from rsfile
       each this_pheno from ind_pheno_cols_ch_gxe_multi
       each chro from list_chro_gemma_gxe
       output:
	 file("${dir_gemma}/${out}.log.txt")
	 set val(our_pheno3), file("${dir_gemma}/${out}.assoc.txt"), val(base) into gemma_manhatten_ch_chro_gxe
       script:
	  our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
	  our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
	  our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
	  data_nomissing     = "pheno-"+our_pheno+".pheno"
	  list_ind_nomissing = "lind-"+our_pheno+".lind"
	  rel_matrix         = "newrel-"+our_pheno+".rel"
	  base               =  plinks[0].baseName
	  inp_fam            =  base+".fam"
	  newbase            =  base+"-"+our_pheno
	  newfam             =  newbase+".fam"
	  gemma_covariate    = "${newbase}.gemma_cov"
	  gemma_gxe          = "${newbase}.gemma_gxe"
	  phef               = "${newbase}_n.phe"
	  covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
	  gxe_opt_gemma      =  (params.gemma_gxe) ? " -gxe $gemma_gxe " : ""
	  out                = "$base-$our_pheno-$chro"
	  dir_gemma          =  (params.gemma_gxe) ? "gemma_gxe" : "gemma"
	  rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
	  """
	  list_ind_nomissing.py --data $covariates --inp_fam $inp_fam --cov_list ${params.covariates},${params.gxe} --pheno $our_pheno3 --dataout $data_nomissing \
				--lindout $list_ind_nomissing
	  gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
	  plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase ${rs_plk_gem}  --chr $chro
	  all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
			     --pheno $our_pheno2 --phe_out ${phef} --form_out 1 --gxe_out $gemma_gxe $gxe_option
	  export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
	  ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf ${params.cut_maf} $gxe_opt_gemma
	  mv output ${dir_gemma}
	  rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
	  """
     }


     gemma_manhatten_ch_chro_gxe_merge=gemma_manhatten_ch_chro_gxe.groupTuple()

     process doMergeGemmaGxE{
          input :
            set (val(this_pheno),file(list_file), base_list) from  gemma_manhatten_ch_chro_gxe_merge
         publishDir "${params.output_dir}/gemma", overwrite:true, mode:'copy'
         output :
             set val(base), val(this_pheno), file("$out") into (gemma_manhatten_ch_gxe_i, gemma_manhatten_ch_gxe)
         script :
             base=base_list[0]
             our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
             our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
             newbase=base+our_pheno
             out = "$base-${our_pheno}.gemma"
             fnames = list_file.join(" ")
             file1  = list_file[0]
             """
             head -1 $file1 > $out
             cat $fnames | grep -v "p_wald" >> $out
             """
     }
   } else {

   ind_pheno_cols_ch = newNamePheno(params.pheno)

   process doGemmaGxE{
    maxForks params.max_forks
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(covariates) from data_ch_gxe
      file(rel) from rel_mat_ch_gxe
      file(plinks) from  gem_ch_gemma_gxe
      file(rsfilelist) from rsfile
    each this_pheno from ind_pheno_cols_ch
    publishDir params.output_dir, overwrite:true, mode:'copy'
    output: 
      file("${dir_gemma}/${out}.log.txt")
      set val(newbase), val(this_pheno), file("${dir_gemma}/${out}.assoc.txt") into (gemma_manhatten_ch_gxe_i, gemma_manhatten_ch_gxe)
    script:
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
       our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
       data_nomissing     = "pheno-"+our_pheno+".pheno" 
       list_ind_nomissing = "lind-"+our_pheno+".lind"
       rel_matrix         = "newrel-"+our_pheno+".rel"
       base               =  plinks[0].baseName
       inp_fam            =  base+".fam"
       newbase            =  base+"-"+our_pheno
       newfam             =  newbase+".fam"
       gemma_covariate    = "${newbase}.gemma_cov"
       gemma_gxe          = "${newbase}.gemma_gxe"
       phef               = "${newbase}_n.phe"
       covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
       gxe_opt_gemma      =  (params.gemma_gxe) ? " -gxe $gemma_gxe " : ""
       out                = "$base-$our_pheno"
       dir_gemma          =  (params.gemma_gxe) ? "gemma_gxe" : "gemma"
       rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
       """
       list_ind_nomissing.py --data $covariates --inp_fam $inp_fam --cov_list ${params.covariates},${params.gxe} --pheno $our_pheno3 --dataout $data_nomissing \
                             --lindout $list_ind_nomissing
       gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase ${rs_plk_gem}
       all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
                          --pheno $our_pheno2 --phe_out ${phef} --form_out 1 --gxe_out $gemma_gxe $gxe_option
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores} 
       ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf ${params.cut_maf} $gxe_opt_gemma
       mv output ${dir_gemma}
       rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
       """
  } 
   }
   
  gemma_manhatten_ch_gxe_freq= gemma_manhatten_ch_gxe_i.combine(Channel.fromPath(params.data, checkIfExists:true)).combine(assoc_ch_gxe_freq)
  process AddedFreqGxEGemma{
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      set val(newbase), val(our_pheno), file(gemmares), file(data_file),file(bed),file(bim),file(fam) from gemma_manhatten_ch_gxe_freq
    publishDir "${params.output_dir}/gemma_gxe", overwrite:true, mode:'copy'
    output:
        file(gemmaresfreq) 
    script :
       base = bed.baseName
       our_pheno2         = our_pheno.replaceAll(/^[0-9]+@@@/,"")
       our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       gemmaresfreq=gemmares+'.withfreq' 
       """
       added_freq_gxe.py --bfile $base --file_gxe $gemmares --pheno_file $data_file --pheno ${our_pheno3} --pheno_gxe ${params.gxe} --out $gemmaresfreq  --plk_cores ${params.gemma_num_cores}
       """
  } 

  process showGemmaManhattanGxE { 
    memory params.other_process_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      set val(base), val(this_pheno), file(assoc) from gemma_manhatten_ch_gxe
    output:
      file("${out}*")  into report_gemma_ch_GxE
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C056$our_pheno"
      """
      general_man.py  --inp $assoc --phenoname $this_pheno --out ${out} --chro_header chr --pos_header ps --rs_header rs --pval_header p_wald --beta_header beta --info_prog "Gemma,GxE: ${params.gxe}"
      """
  }

    
  } else {
    report_gemma_ch_GxE=Channel.empty()
  } 


if (params.assoc+params.fisher+params.logistic+params.linear > 0) {

   process computeTest {
      // template 
      cpus num_assoc_cores
      time params.big_time
      input:
       set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from assoc_ch    
       file (phenof) from pheno_ch
      each test from requested_tests
      each pheno_name from pheno_label_ch
      publishDir "${params.output_dir}/${test}", overwrite:true, mode:'copy'
      output:
        set val(test), val(pheno_name), file("${outfname}.*") into out_ch
      script:
       base = "cleaned"
       pheno_name = pheno_name.replaceFirst("/.*","")
       perm = (params.mperm == 0 ? "" : "--mperm=${params.mperm}")
       adjust = (params.adjust ? "--adjust" : "")
       outfname = "${pheno_name}"
       //test = test_choice 
       if (params.data == "") {
           pheno_cmd = ""
           out = base
       } else {
           pheno_cmd = "--pheno $phenof --pheno-name $pheno_name "
           if (params.covariates) covariate = "--covar ${phenof} --covar-name ${params.covariates} "
           out = pheno
       }
       if(params.sexinfo_available==1){
         sexallow=""
       }else{
        sexallow="--allow-no-sex"
       }
       template "test.sh"
   }


 // log_out_ch = Channel.create()
 
  //log_out_ch.subscribe { println "Completed plink test ${it[0]}" }
 
  process drawPlinkResults { 
    memory params.other_process_mem_req
    input:
    set val(test), val(pheno_name), file(results) from out_ch//.tap(log_out_ch)
    output:
      set file("${base}*man*png"), file ("${base}*qq*png"), file("C050*tex") into report_plink
    publishDir params.output_dir, overwrite:true, mode:'copy'
    script:
      base="cleaned-${test}"
      """
      plinkDraw.py  C050 $base $test ${pheno_name} $gotcovar png
      """
  }

  report_plink_ch=report_plink.groupTuple()

} else {
  report_plink_ch = Channel.empty()
}




if (params.plink_gxe==1) {
  data_ch_plk_gxe = Channel.fromPath(params.data, checkIfExists:true)
  pheno_label_ch_gxe = Channel.from(params.pheno.split(","))
  
  if (params.rs_list=="")
    rsfile_plkgxe=Channel.fromPath("${dummy_dir}/08", checkIfExists:true)
  else
    rsfile_plkgxe=file(params.rs_list, checkIfExists=true)

  process computePlinkGxE {
    cpus num_assoc_cores
    memory plink_mem_req
    time params.big_time
    input:
       set file(filebed),file(filebim),file(filefam) from assoc_ch_gxe
       file (phenof) from data_ch_plk_gxe
       file(rsfile) from rsfile_plkgxe
    each pheno_name from pheno_label_ch_gxe
    publishDir "${params.output_dir}/plink_gxe", overwrite:true, mode:'copy'
    output:
       set file("${out}.qassoc.gxe"),file("${outftmp}.notfind")
       set val(base),val(pheno_name), file("$outf")  into res_plink_gxe
    script:
       pheno_name = pheno_name.replaceFirst("/.*","")
       base       = filebed.baseName
       out        = "$base-${pheno_name}"
       outftmp       = "${out}.tmp.final.gxe"
       outf       = "${out}.qassoc.final.gxe"
       rs_plk        =  (params.rs_list) ?  " --extract  $rsfile" : ""
       """
       PosCol=`head -1 $phenof|sed 's/[\\t ]/\\n/g'|grep -n $params.gxe|awk -F':' '{print \$1-2}'`
       plink --bfile $base --pheno $phenof --pheno-name $pheno_name --threads $num_assoc_cores --out $out --gxe \$PosCol --covar $phenof $rs_plk  --keep-allele-order 
       merge_bim_gxeplink.py --plgxe ${out}.qassoc.gxe --bim $filebim --out $outftmp
       added_freq_gxe.py --bfile $base --file_gxe $outftmp --pheno_file $phenof --pheno ${pheno_name} --pheno_gxe ${params.gxe} --out $outf --plk_cores ${num_assoc_cores} --gwas_chr CHR --gwas_rs SNP


       """
   }

   process showPlinkManhattanGxE {
    label 'gcta'
    memory params.other_process_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      set val(base), val(this_pheno), file(assoc) from res_plink_gxe
    output:
      file("${out}*")  into report_plink_gxe
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C057$our_pheno"
      """
      general_man.py  --inp $assoc --phenoname $this_pheno --out ${out} --chro_header CHR --pos_header POS --rs_header SNP --pval_header P_GXE --beta_header Z_GXE --info_prog "Plink,GxE : ${params.gxe}"
      """
  }
} else {
  report_plink_gxe=Channel.empty()
}

if(params.bgen!="" && (params.fastgwa+params.saige) >0){
     bgen_ch_fastgwa_i=Channel.fromPath(params.bgen, checkIfExists:true)
     process indexbgen {
      label 'utils'
      input :
        path(bgen) from bgen_ch_fastgwa_i
      output :
        path(bgen)  into bgen_ch_fastgwa
        tuple path(bgen), path("${bgen}.bgi") into bgen_ch_saige
        path("${bgen}.bgi") into bgen_ch_fastgwa_idx
      """
      bgenix -g $bgen -index
      """
     }
} else if(params.list_bgen!="" && (params.fastgwa+params.saige) >0){
     bgen_ch_fastgwa_i=Channel.fromPath(file(params.list_bgen).readLines(), checkIfExists:true)
     process indexbgen_list {
      label 'utils'
      input :
        path(bgen) from bgen_ch_fastgwa_i
      output :
        path("${bgen}.bgi") into bgen_ch_idx_1
        tuple path(bgen),path("${bgen}.bgi")  into bgen_ch_saige
      """
      bgenix -g $bgen -index
      """
     }
      
    bgen_ch_fastgwa=Channel.fromPath(file(params.list_bgen).readLines(), checkIfExists:true).collect()
    bgen_ch_fastgwa_idx=bgen_ch_idx_1.collect()
      
}else{
    bgen_ch_fastgwa=Channel.fromPath("${dummy_dir}/06", checkIfExists:true)
    bgen_ch_saige=Channel.fromPath("${dummy_dir}/06", checkIfExists:true).combine(Channel.fromPath("${dummy_dir}/07", checkIfExists:true))
    bgen_ch_fastgwa_idx=Channel.fromPath("${dummy_dir}/08", checkIfExists:true)
}





if(params.fastgwa==1){
   if(params.gcta_grmfile==""){
    process FastGWADoGRM{
       cpus params.fastgwa_num_cores
       maxForks params.max_forks
       label 'gcta'
       memory params.fastgwa_mem_req
       input :
	set file(bed),file(bim),file(fam) from grlm_assoc_ch
	file file_rs from filers_matrel_mat_GWA
      each mpart from 1..params.grm_nbpart
      output : 
       file("mgrm.part_*.grm.id") into idgrm
       file("mgrm.part_*.grm.bin") into bingrm
       file("mgrm.part_*.grm.N.bin") into nbingrm
      script :
	rs_list = balise_filers_rel==1 ? " --extract  $file_rs " : ""

	base   = bed.baseName
	"""
	hostname
	${params.gcta_bin} --bfile $base --make-grm-part  ${params.grm_nbpart} $mpart --thread-num ${params.fastgwa_num_cores} --out mgrm $rs_list --maf  ${params.grm_maf}
	"""
    }


     idgrm_c=idgrm.collect()
     bingrm_c=bingrm.collect()
     nbingrm_c=nbingrm.collect()
     process MergFastGWADoGRM{
       label 'gcta'
       memory params.fastgwa_mem_req
       cpus params.fastgwa_num_cores
       input :
	 file(idgrmallf) from idgrm_c
	 file(bingrmallf) from bingrm_c
	 file(nbingrmallf) from nbingrm_c
       publishDir "${params.output_dir}/fastgwa/grm/", overwrite:true, mode:'copy'
       output :
	 set val(head),file('test_sp_grm.grm.id'), file('test_sp_grm.grm.sp') into grm_all
       script :
	 head='test_sp_grm' 
	 """
	 cat mgrm.part_*_*.grm.id > test_grm.grm.id
	 cat mgrm.part_*_*.grm.bin > test_grm.grm.bin
	 cat mgrm.part_*_*.grm.N.bin > test_grm.grm.N.bin
	 ${params.gcta_bin} --grm test_grm --make-bK-sparse ${params.grm_cutoff} --out $head --thread-num ${params.fastgwa_num_cores}
	 """
     }
   }else{
     grm_all=Channel.from("${params.gcta_grmfile}").combine(Channel.fromPath("${params.gcta_grmfile}.grm.id", checkIfExists:true).combine(Channel.fromPath("${params.gcta_grmfile}.grm.sp", checkIfExists:true)))
   }
   data_ch_fastgwa= Channel.fromPath(params.data, checkIfExists:true)


   pheno_spl_gcta=params.pheno.split(',')

    if (params.covariates)
      covariate_option = "--cov_list ${params.covariates}"
   else
     covariate_option = ""
   covariates_type=""+params.covariates_type
   balqualcov=covariates_type!="" & covariates_type.split(',').contains('1') 
   balquantcov=covariates_type!="" & covariates_type.split(',').contains('0') 


   process FastGWARun{
       maxForks params.max_forks
       label 'gcta'
       memory params.fastgwa_mem_req
       cpus params.fastgwa_num_cores
       input :
	  tuple val(head),path(alldigrm), path(allbingrm) from grm_all
	  tuple path(bed),path(bim),path(fam) from fastgwa_assoc_ch 
	  path(covariates) from data_ch_fastgwa
          path(bgen)  from bgen_ch_fastgwa
          path(bgenindex) from bgen_ch_fastgwa_idx
          path(bgensample) from bgensample_ch_fastgwa
       publishDir "${params.output_dir}/fastgwa/", overwrite:true, mode:'copy'
       each this_pheno from pheno_spl_gcta
       output :
	  tuple val(base), val(this_pheno), path("${out}.fastGWA") into fastgwa_manhatten_ch
       script :
	base=bed.baseName
	phef = "${base}_fastgwa_n.phe"
	covfilequant= "${base}_fastgwa_n.cov"
	covfilequal = "${base}_fastgwa_n.covqual"
	covquant_fastgwa = (balquantcov) ?  " --qcovar $covfilequant " : ""
	our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
	our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
	our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
	out                = "$base-$our_pheno"
	covqual_fastgwa = (balqualcov) ? " --covar $covfilequal " : ""
	covqual_cov = (balqualcov) ? " --cov_type ${params.covariates_type} --covqual_file $covfilequal " : ""
        genet=(params.bgen=='') ? "  --bfile $base " : " --bgen $bgen --sample $bgensample --info ${params.bgen_mininfo} "
        listbgen="list_bgen"
        genet=(params.list_bgen=='') ? "  $genet " : " --mbgen $listbgen --sample $bgensample --info ${params.bgen_mininfo} "
        baliselbgen=(params.list_bgen=='') ? 0 : 1
          
	"""
        if [ "$baliselbgen" -eq "1" ]
        then
        ls *.bgen > $listbgen
        fi
	all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option --pheno ${this_pheno} --phe_out ${phef}  --cov_out $covfilequant --form_out 4  $covqual_cov --nona 1
	${params.gcta_bin} $genet ${params.fastgwa_type}  --pheno $phef  $covquant_fastgwa --threads ${params.fastgwa_num_cores} --out $out --grm-sparse $head $covqual_fastgwa --maf ${params.cut_maf}
	"""
   }

   process showFastGWAManhattan {
   label 'py3fast'
   memory params.other_process_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      set val(base), val(this_pheno), file(assoc) from fastgwa_manhatten_ch
    output:
      file("${out}*")  into report_fastgwa_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C052-fastGWA-"+our_pheno
      //CHR	SNP	POS	A1	A2	N	AF1	BETA	SE	P
      """
      general_man.py  --inp $assoc --phenoname $this_pheno --out ${out} --chro_header CHR --pos_header POS --rs_header SNP --pval_header P --beta_header BETA --info_prog fastGWA
      """
  }

}else{
   report_fastgwa_ch=Channel.empty()
}

if(params.saige==1){
    bim_ch_saige=channel.fromPath(bim)
    process getListeChro_saige{
        input :
          file(BimFile) from bim_ch_saige
        output :
          stdout into chrolist_saige
        script:
         """
         cat $BimFile|awk '{print \$1}'|uniq|sort|uniq
        """
   }
   check3 = Channel.create()
   chrolist_saige_ch=chrolist_saige.flatMap { list_str -> list_str.split() }.tap ( check3)


  if(balise_filers_rel==1){
    process subplink_heritability_saige{
    input :
       file(rs) from filers_her_saige
       file(plk) from ch_saige_heritability
    output :
       set file("${out}.bed"), file("${out}.bim"), file("${out}.fam")  into ch_plk_saige
    script :
       bfile=plk[1].baseName
       out=bfile+'_subrs'
 //--extract $rs
       """
       plink -bfile $bfile --extract $rs --make-bed -out $out --keep-allele-order ${params.snp_rel_param_plk} --allow-extra-chr 
       """
    }
   }else{
   ch_plk_saige=Channel.create()
    Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { ch_plk_saige}

  }


  fam_ch_saige = Channel.create()
  plk_ch_saige_her = Channel.create()
  ch_plk_saige.separate (plk_ch_saige_her, fam_ch_saige) { a -> [ a, a[2]] }
    data_ch_saige = Channel.fromPath(params.data,checkIfExists:true)

  
  if(params.list_vcf!=''){
   firstvcf_ch=Channel.fromPath(file(params.list_vcf).readLines()[0], checkIfExists:true)
   process checkidd_saige_vcf{
    label 'R'
    input :
      path(covariates) from data_ch_saige
      path(vcf) from  firstvcf_ch
      tuple path(bed), path(bim), path(fam) from plk_ch_saige_her
    output :
      tuple path("${bfileupdate}.bed"), path("${bfileupdate}.bim"), path("${bfileupdate}.fam") into plk_ch_saige_her_idupd
      tuple path(covariates_form),path("${bfileupdate}.fam") into data_ch_saige_form
    script :
      covariates_form=covariates.baseName+'_saige.ind'
      bfile=bed.baseName
      bfileupdate=bfile+'_idsaige'
      """
      zcat $vcf |head -1000 |grep "#"| tail -1|awk '{for(Cmt=10;Cmt<=NF;Cmt++)print \$Cmt}' > fileind
      awk '{print \$1"\t"\$1}' fileind > keep
      format_saige_pheno.r --data $covariates --ind_vcf fileind --out ${covariates_form} --pheno ${params.pheno}  --pheno_bin ${params.pheno_bin}
      awk '{if( \$1 ~ /^[0-9+]\$/)print \$1"\t"\$4"\t"\$4"\t"\$4}' $bim > keep.range
      plink -bfile  $bfile --keep-allele-order -out $bfileupdate --make-bed --update-ids  ${covariates_form}"_updateid" --keep keep --allow-extra-chr --extract range keep.range
      """
   }
  }else{
   process checkidd_saige{
    label 'R'
    input :
      path(covariates) from data_ch_saige
      path(bgensample) from  bgensample_ch2
      tuple path(bed), path(bim), path(fam) from plk_ch_saige_her
    output :
      tuple path("${bfileupdate}.bed"), path("${bfileupdate}.bim"), path("${bfileupdate}.fam") into plk_ch_saige_her_idupd
      tuple path(covariates_form),path("${bfileupdate}.fam") into data_ch_saige_form
    script :
      covariates_form=covariates.baseName+'_saige.ind'
      bfile=bed.baseName
      bfileupdate=bfile+'_idsaige'
      indnbgen=(params.bgen!='') ? "--ind_bgen $bgensample" : ""
      """
      format_saige_pheno.r --data $covariates $indnbgen --out ${covariates_form}  --pheno ${params.pheno}  --pheno_bin ${params.pheno_bin}
      awk '{if( \$0 ~ /^[0-9+]\$/)print \$1"\t"\$4"\t"\$4"\t"\$4}' $bim > keep.range
      plink -bfile  $bfile --keep-allele-order -out $bfileupdate --make-bed --update-ids  ${covariates_form}"_updateid" --extract range keep.range
      """
   }
 }
  process  getSaigePheno {
    input:
      tuple path(covariates), path(fam) from data_ch_saige_form
    output:
      file(phef) into saige_data_format_ch
      stdout into pheno_cols_ch_saige
    script:
      covoption = (params.covariates=="") ? "" : " --cov_list ${params.covariates}"
      base = fam.baseName
      phef = "${base}_saige_n.phe"
      """
      all_covariate.py --data  $covariates --inp_fam  $fam $covoption \
                          --pheno ${params.pheno} --phe_out ${phef}  --form_out 2 --nona 1
      """
  }
  covariable_saige=regenieCofact(params.covariates, params.covariates_type,0)

  pheno_ch_saige = Channel.create()
  check = Channel.create()
  pheno_cols_ch_saige.flatMap { list_str -> list_str.split() }.tap ( check) .set { pheno_ch_saige}

  process saige_computed_variance{ 
   memory params.saige_mem_req
   cpus params.saige_num_cores
   label 'saige'
   input :
      tuple path(bed), path(bim), path(fam) from plk_ch_saige_her_idupd
      path(pheno) from saige_data_format_ch
   publishDir "${params.output_dir}/saige/varexp/", overwrite:true, mode:'copy'
   each this_pheno from pheno_ch_saige
   output :
     set val(our_pheno) ,file("${output}.rda"), file("${output}.varianceRatio.txt"), val(plk) into file_varexp 
     path("run_step_1.sh") 
   script :
     //covoption = (params.covariates=="") ? "" : " --covarColList=${params.covariates}"
     binpheno = (params.pheno_bin==1) ? " --traitType=binary " : " --traitType=quantitative"
     Loco = (params.saige_loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     plk=bed.baseName
     our_pheno    = this_pheno.replaceAll(/|\/np.\w+/,"").replaceAll(/[0-9]+@@@/,"")
     output=our_pheno+"_var"
     covoption= (params.covariates=="") ? "" : "  ${covariable_saige} "
     """
     ${params.saige_bin_fitmodel} \
        --plinkFile=$plk \
        --phenoFile=$pheno \
        --phenoCol=$our_pheno $covoption \
        --sampleIDColinphenoFile=IID \
        --outputPrefix=./$output \
        --nThreads=${params.saige_num_cores} $Loco \
        --minMAFforGRM=${params.grm_maf} $binpheno \
        ${params.saige_otheropt_step1} ${params.saige_otheropt}  --maxMissing ${params.cut_geno}
     cp .command.sh run_step_1.sh
     """
  }
  if(params.list_vcf!=''){
   listvcf_ch=Channel.fromPath(file(params.list_vcf).readLines(), checkIfExists:true)
   process buildindex {
   memory { strmem(params.other_process_mem_req) + 5.GB * (task.attempt -1) }               
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }         
   maxRetries 3                                                           
    label 'utils'
    input :
       file(vcf) from listvcf_ch
    output :
        set file(vcf), file(vcfindex) into listvcf_ind_ch
    script :
       vcfindex=vcf+".csi"
       """
       hostname
       tabix -C -p vcf $vcf
       """
   }
   vcf_andvariance_ch=listvcf_ind_ch.combine(file_varexp) 
   process doSaige{
    memory params.saige_mem_req
    cpus params.saige_num_cores
    label 'saige'
    input :
       set file(vcf), file(vcfindex), val(our_pheno) , file(rda), file(varRatio), val(base) from vcf_andvariance_ch 
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.sh"
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.log"
    output :
      set val(our_pheno),file("$output"), val(base) into ch_saige_bychro
      tuple path("${output}.sh"), path("${output}.log")
    script :
     output=vcf.baseName+".res"
     //bin_option_saige= (params.pheno_bin==1) ? " --IsOutputAFinCaseCtrl=TRUE  --IsOutputNinCaseCtrl=TRUE --IsOutputHetHomCountsinCaseCtrl=TRUE " : ""
     Loco = (params.saige_loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     imputed=(params.saige_imputed_data==1) ? "--is_imputed_data=TRUE " : '--is_imputed_data=FALSE' 
     moredetail=(params.pheno_bin==1) ? " --is_output_moreDetails=TRUE " : ""
     """
      Chro=`zcat $vcf|grep -v "#"|head -1|awk '{print \$1}'`
      ${params.saige_bin_spatest} \
        --vcfFile=$vcf\
        --vcfFileIndex=$vcfindex \
        --vcfField=${params.vcf_field} \
        --minMAF=${params.cut_maf}\
        --minMAC=${params.vcf_minmac} \
        --chrom=\$Chro \
        --GMMATmodelFile=$rda \
        --varianceRatioFile=$varRatio \
        --SAIGEOutputFile=$output $Loco $imputed  $moredetail ${params.saige_otheropt_step2}  ${params.saige_otheropt}    --maxMissing ${params.cut_geno}
     cp .command.sh $output".sh"
     cp .command.out $output".log"
     """
   }
 }else if(params.bgen!=''){
   bgen_andvariance_ch=bgen_ch_saige.combine(bgensample_ch_saige).combine(file_varexp)
   process doSaigeBgen{
    memory params.saige_mem_req
    cpus params.saige_num_cores
    label 'saige'
    input :
       tuple path(bgen), path(bgenindex), path(bgensample),val(our_pheno),path(rda), path(varRatio), val(base) from bgen_andvariance_ch
    each Chro from chrolist_saige_ch
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.sh"
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.log"
    output :
      set val(our_pheno),file("$output"), val(base) into ch_saige_bychro
      tuple path("${output}.sh"), path("${output}.log")
    script :
     output=bgen.baseName+"_"+Chro+".saige"
     //bin_option_saige= (params.pheno_bin==1) ? " --IsOutputAFinCaseCtrl=TRUE  --IsOutputNinCaseCtrl=TRUE --IsOutputHetHomCountsinCaseCtrl=TRUE " : ""
     Chro=Chro.replace(' ','')
     Chro=(Chro=="23") ? "X": "$Chro"
     Loco = (params.saige_loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     imputed=(params.saige_imputed_data==1) ? "--is_imputed_data=TRUE --impute_method=${params.saige_impute_method} " : '--is_imputed_data=FALSE' 
     moredetail=(params.pheno_bin==1) ? " --is_output_moreDetails=TRUE " : ""
     """
      ${params.saige_bin_spatest} \
        --bgenFile=$bgen    \
        --bgenFileIndex=$bgenindex \
        --sampleFile=$bgensample \
        --AlleleOrder=ref-first \
        --chrom=$Chro \
        --minMAF=${params.cut_maf}\
        --minMAC=${params.vcf_minmac} \
        --GMMATmodelFile=$rda \
	--SAIGEOutputFile=$output $imputed $Loco $moredetail  ${params.saige_otheropt_step2}   ${params.saige_otheropt}  --maxMissing ${params.cut_geno}
        cp .command.sh run_step_2.sh
     cp .command.out $output".log"
     """
   }

 } else if(params.list_bgen!=''){
   process getchrobgen{
     memory params.saige_mem_req
      label 'utils'
      input :
       tuple path(bgen), path(bgenindex) from bgen_ch_saige
      output :
       tuple env(Chro), path(bgen), path(bgenindex) into bgen_ch_saige_chro
       """
       Chro=`bgenix -g $bgen -list |sed '1,2d'| awk '{print \$3}' |head -1`
       """

   }
   bgen_andvariance_ch=bgen_ch_saige_chro.combine(bgensample_ch_saige).combine(file_varexp)
   process doSaigeListBgen{
    memory params.saige_mem_req
    cpus params.saige_num_cores
    label 'saige'
    input :
       tuple val(Chro),path(bgen), path(bgenindex), path(bgensample),val(our_pheno),path(rda), path(varRatio), val(base) from bgen_andvariance_ch
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.sh"
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.log"
    output :
      set val(our_pheno),file("$output"), val(base) into ch_saige_bychro
      tuple path("${output}.sh"), path("${output}.log")
    script :
     output=bgen.baseName+"_"+Chro+".saige"
     //bin_option_saige= (params.pheno_bin==1) ? " --IsOutputAFinCaseCtrl=TRUE  --IsOutputNinCaseCtrl=TRUE --IsOutputHetHomCountsinCaseCtrl=TRUE " : ""
     Loco = (params.saige_loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     imputed=(params.saige_imputed_data==1) ? "--is_imputed_data=TRUE --impute_method=${params.saige_impute_method}" : '--is_imputed_data=FALSE' 
     moredetail=(params.pheno_bin==1) ? " --is_output_moreDetails=TRUE " : ""
     """
      ${params.saige_bin_spatest} \
        --bgenFile=$bgen    \
        --bgenFileIndex=$bgenindex \
        --sampleFile=$bgensample \
        --chrom=$Chro \
        --AlleleOrder=ref-first \
        --minMAF=${params.cut_maf}\
        --minMAC=${params.vcf_minmac} \
        --GMMATmodelFile=$rda \
        --varianceRatioFile=$varRatio \
        --SAIGEOutputFile=$output ${Loco} $imputed  $moredetail ${params.saige_otheropt_step2} ${params.saige_otheropt}
        cp .command.sh run_step_2.sh
     cp .command.out $output".log"
     """
   }


   

 }else{

   plink_andvariance_ch=ch_saige_assoc.combine(file_varexp)
   process doSaigePlink{
    memory params.saige_mem_req
    cpus params.saige_num_cores
    label 'saige'
    input :
       tuple file(bed), file(bim), file(fam),val(our_pheno) , path(rda), path(varRatio), val(base) from plink_andvariance_ch
    each Chro from chrolist_saige_ch
    output :
      tuple val(our_pheno),file("$output"), val(base) into ch_saige_bychro
      tuple path("${output}.sh"), path("${output}.log")
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.sh"
   publishDir "${params.output_dir}/saige/log/", overwrite:false, mode:'copy', pattern: "*.log"
    script :
     output=bed.baseName+'_'+Chro+".saige"
     Loco = (params.saige_loco==1) ? " --LOCO=TRUE " : " --LOCO=FALSE "
     //bin_option_saige= (params.pheno_bin==1) ? " --IsOutputAFinCaseCtrl=TRUE  --IsOutputNinCaseCtrl=TRUE --IsOutputHetHomCountsinCaseCtrl=TRUE " : ""
     imputed=(params.saige_imputed_data==1) ? "--is_imputed_data=TRUE --impute_method=${params.saige_impute_method} " : '--is_imputed_data=FALSE' 
     moredetail=(params.pheno_bin==1) ? " --is_output_moreDetails=TRUE " : ""
     """

      ${params.saige_bin_spatest} \
        --bedFile=$bed \
        --bimFile=$bim \
        --famFile=$fam \
        --chrom=$Chro \
        --AlleleOrder=alt-first \
        --minMAF=${params.cut_maf}\
        --minMAC=${params.vcf_minmac} \
        --GMMATmodelFile=$rda \
        --varianceRatioFile=$varRatio \
        --SAIGEOutputFile=$output $Loco $imputed  $moredetail ${params.saige_otheropt_step2}  ${params.saige_otheropt}
        cp .command.sh run_step_2.sh
     cp .command.out $output".log"
     """
   }


 }
 ch_saige_res=ch_saige_bychro.groupTuple()
 process doMergeSaige{
          input :
            tuple (val(this_pheno),file(list_file), base_list) from ch_saige_res
         publishDir "${params.output_dir}/saige", overwrite:true, mode:'copy'
         output :
             tuple val(base), val(our_pheno2), file("$out") into saige_manhatten_ch
         script :
             base=base_list[0]
             our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
             our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
             out = "$base-${our_pheno}.saige"
             fnames = list_file.join(" ")
             file1  = list_file[0]
             """
             head -1 $file1 > $out
             cat $fnames | grep -v CHR >> $out
             """
  }
  process showSaigeManhattan {
    memory params.other_process_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      set val(base), val(this_pheno), file(assoc) from saige_manhatten_ch
    output:
      file("${out}*")  into report_saige_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C058-saige-"+our_pheno
      """
      general_man.py  --inp $assoc --phenoname $this_pheno --out ${out} --chro_header CHR --pos_header POS --rs_header  SNPID,MarkerID --pval_header p.value --beta_header BETA --info_prog SAIGE
      """
  }



}else{
 report_saige_ch=Channel.empty()
}

   def regenieCofact(args,infoargs, regenie) {
      //Method code
      infoargs=""+infoargs
      splargs=args.split(",")
      splinfoargs=infoargs.split(",")
      if(splinfoargs.size()==0){
        cov="--covarColList "+args
        return(cov)
      }
      if(splargs.size() != splinfoargs.size()){
         System.err.println("args and args type for Boltlmm was not same size : "+args+" "+infoargs)
         System.exit(-11)
      }
      Cofactqual=""
      Cofactquant=""
      allcov=""
      for (i = 0; i <splargs.size(); i++) {
          /*0 : for quantitatif */
          /* 1 for qualitatif*/
          if(allcov=="")allcov=splargs[i]
          else allcov+=","+splargs[i]
          if     (splinfoargs[i]=='1'){
              if(Cofactqual=="")Cofactqual=splargs[i]
              else Cofactqual+=","+splargs[i]

          }else {
           if(splinfoargs[i]=='0'){
              if(Cofactquant=="")Cofactquant=splargs[i]
              else Cofactquant+=","+splargs[i]
           }else{
             System.err.println("type args for "+splargs[i]+" doesn't know "+ splinfoargs[i]+"\n 1 for quantitatif arguments\n 0 for qualitatif arguments")
             System.exit(-10)
          }
        }
      }
     CofactStr=""
     if(regenie==1){
       if(Cofactqual!="")CofactStr=" --catCovarList "+Cofactqual
        if(Cofactquant!="") CofactStr+=" --covarColList "+allcov
     }else{
       if(Cofactqual!="")CofactStr=" --qCovarColList  "+Cofactqual
        CofactStr+=" --covarColList "+allcov
     }
      return(CofactStr)
   }


                                                                        

if(params.regenie==1){

 data_ch_regenie = Channel.fromPath(params.data, checkIfExists:true)
 check = Channel.create()
 pheno_cols_ch_regenie= Channel.from(params.pheno.split(","))

  if(params.bgen!=""){
    bgen_ch_regenie=Channel.fromPath(params.bgen, checkIfExists:true)
    if(params.bgen_sample==''){
      println "params.bgen_sample not initialise when bgen params initial";
      System.exit(-2);
     }
  }else if(params.list_bgen!="") {
    bgen_ch_regenie=Channel.fromPath(file(params.list_bgen).readLines(), checkIfExists:true)
    if(params.bgen_sample==''){
      println "params.bgen_sample not initialise when bgen params initial";
      System.exit(-2);
     }

   }else{
    bgen_ch_regenie=Channel.fromPath("${dummy_dir}/06", checkIfExists:true)
  }




 covariable_regenie=regenieCofact(params.covariates, params.covariates_type,1)
   

 process regenie_step1{
   cpus params.regenie_num_cores
   memory params.regenie_mem_req
   input :
     path(data) from data_ch_regenie
     path(rsrel) from filers_matrel_regenie
     tuple path(bed), path(bim), path(fam) from ch_regenie_assoc
   each pheno from pheno_cols_ch_regenie
   //errorStrategy { task.exitStatus == 1 ; return 'ignore' }
   publishDir "${params.output_dir}/regenie/step1", overwrite:true, mode:'copy', pattern: "*.loco"
   publishDir "${params.output_dir}/regenie/step1", overwrite:true, mode:'copy', pattern: "*.list"
   publishDir "${params.output_dir}/regenie/step1", overwrite:true, mode:'copy', pattern: "*.report"
   output : 
    tuple val(our_pheno), path("$phef"),path("${out}_pred.list"), path("${out}_1.loco"),path(bed), path(bim), path(fam),  optional :true into ch_regenie_pheno
    path("*.report")
    path("*.log")
   script :
       our_pheno       = pheno.replaceAll(/\/np.\w+/,"").replaceAll(/[0-9]+@@@/,"")
      phef=pheno+".pheno"
      loco=(params.regenie_loco==0) ? "" : " --loocv "
      bfile=bed.baseName
      bfilesub=bfile+"_sub"
      keeppos=(rsrel.toString()=='00') ? ""   : " --extract $rsrel "
      print "keeppos $keeppos"
      covoption = (params.covariates=="") ? "" : " --cov_list ${params.covariates}"
      covoption_regenie= (params.covariates=="") ? "" : " --covarFile $phef ${covariable_regenie} "
      bsize=(params.regenie_bsize_step1==0) ? " ${params.regenie_bsize} " : " ${params.regenie_bsize_step1}"
      regenie_loco=(params.regenie_loco=="") ? "" : " --loocv "
      out=phef+"_regenie"
      println covoption_regenie
      gxe=(params.regenie_gxe==0) ? "" : " --gxe ${params.gxe} "
      """
      all_covariate.py --data  $data --inp_fam  $fam  $covoption \
                          --pheno $pheno --phe_out ${phef}  --form_out 2 --nona  1 $gxe
      plink -bfile $bfile $keeppos --make-bed -out $bfilesub --keep $phef ${params.snp_rel_param_plk}
      ${params.regenie_bin} --step 1   --bed $bfilesub    --phenoFile $phef  --phenoCol ${our_pheno} --bsize $bsize $regenie_loco --out  $out --threads ${params.regenie_num_cores} ${params.regenie_otheropt_step1} $covoption_regenie
      if [ ! -f $out"_1.loco" ]
      then
      touch $out"_1.loco"
      fi
      cp .command.sh "${our_pheno}"_regenie_step1.cmd.report
      cp .command.log "${our_pheno}"_regenie_step1.log.report
      cp .command.err "${our_pheno}"_regenie_step1.err.report
      """
 }
 ch_regenie_pheno_2=ch_regenie_pheno.combine(bgen_ch_regenie).combine(bgensample_ch_regenie)
 process regenie_step2{
   cpus params.regenie_num_cores
   memory params.regenie_mem_req
   input :
    tuple val(pheno), path(data),path(list), path(loco),path(bed), path(bim), path(fam),path(bgen), path(bgensample)  from ch_regenie_pheno_2
   publishDir "${params.output_dir}/regenie/step2", overwrite:true, mode:'copy', pattern: "*.cmd"
   output :
     tuple val(pheno), path("${out}*${pheno}.regenie"),val(bfile) into regenie_manhatten_chi
     path("*regenie_step2.cmd")
   script :
    loco=(params.regenie_loco==0) ? "" : " --loocv "
    bfile=bed.baseName
    covoption_regenie= (params.covariates=="") ? "" : " --covarFile $data   ${covariable_regenie} "
    bsize=(params.regenie_bsize_step2==0) ? " ${params.regenie_bsize} " : " ${params.regenie_bsize_step2}"
    regenie_loco=(params.regenie_loco=="") ? "" : " --loocv "
    genet=(params.bgen=="")? " --bed $bfile " : " --bgen $bgen --sample $bgensample "
    genet=(params.list_bgen=="")? " $genet " : " --bgen $bgen --sample $bgensample "
    loco=(params.regenie_loco==0) ? "" : " --loocv "
    out=pheno+"_regenie_assoc"
    out=(params.list_bgen=="")? out: bgen.baseName
    gxe=(params.regenie_gxe==0) ? "" : " --interaction ${params.gxe} "
    """
     ${params.regenie_bin} --step 2  $genet   --phenoFile $data --phenoCol $pheno ${covoption_regenie} --bsize $bsize   --pred $list  $loco   --out $out --threads ${params.regenie_num_cores} ${params.regenie_otheropt_step2} $gxe  
     cp .command.sh "${pheno}"_regenie_step2.cmd
  """
 }

 regenie_manhatten_chi_m=regenie_manhatten_chi.groupTuple()
 process merge_regenie{
    input :
       set (val(pheno),file(list_file), base_list) from  regenie_manhatten_chi_m
    publishDir "${params.output_dir}/regenie/", overwrite:true, mode:'copy'
    output :
      tuple val(bfile), val(pheno), path(assoc) into regenie_manhatten_ch_merge
    script :
      bfile=base_list[0]
      assoc= "${pheno}.regenie"
      fnames = list_file.join(" ")
      file1  = list_file[0]
      """
      head -1 $file1 > $assoc
      cat $fnames | grep -v CHROM >> $assoc
      """
 }

 process format_regeniesumstat{
  input :
   tuple val(bfile), val(pheno), path(assoc) from regenie_manhatten_ch_merge
  publishDir "${params.output_dir}/regenie/", overwrite:true, mode:'copy'
  output :
    tuple val(bfile), val(pheno), path(newassoc) into regenie_manhatten_ch
  script :
    newassoc=assoc.baseName+"_format.regenie" 
  """
  format_sumstat_regenie.py $assoc $newassoc
  """ 
 }

  process showRegenieManhattan {
    memory params.other_process_mem_req
    publishDir params.output_dir, overwrite:true, mode:'copy'
    input:
      set val(base), val(this_pheno), file(assoc) from regenie_manhatten_ch
    output:
      file("${out}*")  into report_regenie_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C059-regenie-"+our_pheno
      //CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
      """
      general_man.py  --inp $assoc --phenoname $this_pheno --out ${out} --chro_header CHROM --pos_header GENPOS --rs_header  ID --pval_header P --beta_header BETA --info_prog regenie
      """
  }



}else{

 report_regenie_ch=Channel.empty()

}




def getres(x) {
  def  command1 = "$x"
  def  command2 = "head -n 1"
  def proc1 = command1.execute()
  def proc2 = command2.execute()
  def proc = proc1 | proc2
  proc.waitFor()              
  res ="${proc.in.text}"
  return res.trim()
}

//nextflowversion =getres("nextflow -v")
nextflowversion =nextflow.version
if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"

report_ch = report_fastlmm_ch.flatten().mix(pheno_report_ch.flatten())
                                     .mix(report_pca_ch.flatten())
				     .mix(report_plink_ch.flatten())
				     .mix(report_bolt_ch.flatten())
				     .mix(report_fastgwa_ch.flatten())
				     .mix(report_gemma_ch_GxE.flatten())
				     .mix(report_plink_gxe.flatten())
				     .mix(report_saige_ch.flatten())
                                     .mix(report_regenie_ch.flatten())
                                     .mix(report_gemma_ch.flatten()).toList()

process doReport {
  label 'latex'
  input:
    file(reports) from report_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file("${out}.pdf")
  script:
    out = params.output+"-report"
    these_phenos     = params.pheno
    these_covariates = params.covariates
    config = getConfig()
    images = workflow.container
    texf   = "${out}.tex"
    template "make_assoc_report.py"
}





