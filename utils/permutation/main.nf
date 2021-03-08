#!/usr/bin/env nextflow
/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Shaun Aron
 *   	Rob Clucas
 *      Eugene de Beste
 *      Lerato Magosi
 *      Jean-Tristan Brandenburg
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths




def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates","gemma_num_cores","gemma_mem_req","gemma","linear","logistic","assoc","fisher", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "other_mem_req", "max_plink_cores", "pheno","big_time","thin", "gemma_mat_rel","print_pca", "file_rs_buildrelat","genetic_map_file", "rs_list","adjust","mperm"]


/*JT : append argume boltlmm, bolt_covariates_type */
/*bolt_use_missing_cov --covarUseMissingIndic : “missing indicator method” (via the --covarUseMissingIndic option), which adds indicator variables demarcating missing status as additional covariates. */
//ParamBolt=["bolt_ld_scores_col", "bolt_ld_score_file","boltlmm", "bolt_covariates_type",  "bolt_use_missing_cov", "bolt_num_cores", "bolt_mem_req", "exclude_snps", "bolt_impute2filelist", "bolt_impute2fidiid", "bolt_otheropt"]
//allowed_params+=ParamBolt
//ParamFast=["fastlmm","fastlmm_num_cores", "fastlmm_mem_req", "fastlmm_multi", "fastlmmc_bin"]
//allowed_params+=ParamFast
/*Gxe : */
GxE_params=['gemma_gxe', "plink_gxe", "gxe"]
allowed_params+=GxE_params


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
params.print_pca = 1
params.file_rs_buildrelat = ""
params.genetic_map_file = ""
outfname = params.output_testing



/* Defines the path where any scripts to be executed can be found.
 */

/* Do permutation testing -- 0 for none, otherwise give number */
params.mperm = 1000

/* Adjust for multiple correcttion */
params.adjust = 0

supported_tests = ["assoc","fisher","model","cmh","linear","logistic","boltlmm", "fastlmm", "gemma", "gemma_gxe"]


params.assoc     = 0
params.fisher   = 0
params.cmh     =  0
params.model   =  0
params.linear   = 0
params.logistic = 0
params.gemma = 0
params.gemma_mem_req = "6GB"
params.gemma_relopt = 1
params.gemma_lmmopt = 4
params.gemma_mat_rel = ""
params.gemma_num_cores = 8
params.pheno = "_notgiven_"

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
params.fastlmmc_bin =""
/*gxe param : contains column of gxe*/
params.gemma_gxe=0
params.plink_gxe=0
params.max_plink_cores = 4
params.rs_list=""
params.gxe=""


params.input_pat  = 'raw-GWA-data'

params.sexinfo_available = "false"


params.plink_mem_req = '750MB' // how much plink needs for this
params.other_process_memory = '750MB' // how much other processed need


plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_memory
max_plink_cores = params.max_plink_cores 

params.help = false


data_ch = file(params.data)

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
  f = new File(fname)
  if (! f.exists()) {
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
fileColExists(params.data,"${params.gxe} - gxe", params.gxe)

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




bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()


gemma_assoc_ch = Channel.create()
/*JT initatilisation of boltlmm_assoc_ch*/
boltlmm_assoc_ch = Channel.create()
fastlmm_assoc_ch = Channel.create()
rel_ch_fastlmm = Channel.create()

pca_in_ch = Channel.create()
assoc_ch  = Channel.create()
assoc_ch_gxe  = Channel.create()
raw_src_ch= Channel.create()

Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_src_ch }


println "\nTesting data            : ${params.input_pat}\n"
println "Testing for phenotypes  : ${params.pheno}\n"
println "Using covariates        : ${params.covariates}\n\n"

if (params.gemma) println "Doing gemma testing"
if (params.assoc) println "Doing assoc testing"
if (params.linear) println "Doing linear regression testing"
if (params.logistic) println "Doing logistic regression testing"
if (params.fastlmm == 1) println "Doing mixed model with fastlmm "
if (params.boltlmm == 1) println "Doing mixed model with boltlmm "
if(params.gemma_gxe==1)println "Doing mixed model with gemma and gxe with "+params.gxe
if(params.plink_gxe==1)println "Doing with plink gxe with "+params.gxe
println "\n"

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
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  ( pca_in_ch, assoc_ch, gemma_assoc_ch, boltlmm_assoc_ch,fastlmm_assoc_ch, rel_ch_fastlmm)
    script:
       base = bed.baseName
       out  = base+"_t"
       "plink --keep-allele-order --bfile $base $thin $chrom --make-bed --out $out"
  }

  println "\nData has been thinned or only some chromosomes used  (is the run for test purposes only?)\n"
   


} else {
    /*JT : append boltlmm_assoc_ch and a]*/
    raw_src_ch.separate( pca_in_ch, assoc_ch, assoc_ch_gxe, gemma_assoc_ch, boltlmm_assoc_ch, fastlmm_assoc_ch,rel_ch_fastlmm) { a -> [a,a,a,a,a,a,a] }
}





num_assoc_cores = params.mperm == 0 ? 1 : Math.min(10,params.max_plink_cores)

supported_tests = ["assoc","fisher","model","cmh","linear","logistic"]

requested_tests = supported_tests.findAll { entry -> params.get(entry) }


covariate = ""
gotcovar  = 0
pheno     = ""




  def newNamePheno(Pheno){
      SplP=Pheno.split(',')
      for (i = 0; i <SplP.size(); i++) {
         SplP[i]=(i+1)+"-"+SplP[i]
      }
      return(SplP)
  }




 
if (params.data != "") {

   checker(file(params.data))

   if (params.covariates != "") {
      gotcovar = 1
  }


  
  process extractPheno {
    input:
     file(data) from data_ch
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
    input:
    file(data) from data_ch
    output:
      file ("B050*") into report_ch
    script:
      "phe_distrib.py --pheno ${params.pheno} $data B050 "
  }
}  else {
  report_ch = Channel.empty()
  pheno_label = ""
  pheno_label_ch = Channel.from("")
}



if (params.gemma+params.gemma_gxe>0) {
   if(params.file_rs_buildrelat==""){
        filers_matrel_mat_gem=file('NO_FILE')
     }else{
        filers_matrel_mat_gem=Channel.fromPath(params.file_rs_buildrelat)
   }

  rel_ch_gemma = Channel.create()
  gem_ch_gemma = Channel.create()
  gem_ch_gemma_gxe = Channel.create()
  gemma_assoc_ch.separate (rel_ch_gemma, gem_ch_gemma, gem_ch_gemma_gxe) { a -> [a, a, a] }
  if(params.gemma_mat_rel==""){
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
       gemma -bfile $base  -gk ${params.gemma_relopt} -o $base -p pheno -n 3 $rs_list
       """
  }
  }else{
   rel_mat_ch=Channel.fromPath(params.gemma_mat_rel) 
   rel_mat_ch_gxe=Channel.fromPath(params.gemma_mat_rel) 
  }
}

if (params.gemma == 1){

  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""
  ind_pheno_cols_ch = newNamePheno(params.pheno)
   if(params.rs_list==""){
        rsfile=file('NO_FILE5')
     }else{
        rsfile=file(params.rs_list)
   }

  process doFormatGemma{
    input:
      file(covariates) from data_ch
      file(rel) from rel_mat_ch
      file(plinks) from  gem_ch_gemma
      file(rsfilelist) from rsfile
    each this_pheno from ind_pheno_cols_ch
    output:
      set file("${newbase}.bed"),file("${newbase}.bim"),file("${newbase}.fam") , file(rel_matrix),val(this_pheno), file(gemma_covariate), file(phef) into (gemma_data_perm1, gemma_data_perm2)
    script:
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+-/,"")
       ourpheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
       data_nomissing     = "pheno-"+our_pheno+".pheno"
       list_ind_nomissing = "lind-"+our_pheno+".lind"
       rel_matrix         = "newrel-"+our_pheno+".rel"
       base               =  plinks[0].baseName
       inp_fam            =  base+".fam"
       newbase            =  base+"-"+our_pheno
       newfam             =  newbase+".fam"
       gemma_covariate    = "${newbase}.gemma_cov"
       phef               = "${newbase}_n.phe"
       rs_plk_gem         =  (params.rs_list) ?  " --extract  $rsfilelist" : ""
       """
       list_ind_nomissing.py --data $covariates --inp_fam $inp_fam $covariate_option --pheno $ourpheno3 --dataout $data_nomissing \
                             --lindout $list_ind_nomissing
       gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem}
       all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
                          --pheno $our_pheno2 --phe_out ${phef} --form_out 1
       """



  }

  process doGemmaPerm{
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      set file(bedfile),file(bimfile),file(famfile) , file(rel_matrix),val(this_pheno), file(gemma_covariate), file(phefi) from gemma_data_perm1
    each mperm from 1..params.mperm
    publishDir params.output_dir, overwrite:true, mode:'copy'
    output:
      set val(our_pheno), file("${dir_gemma}/${this_pheno}/${out}.assoc.txt") into gemma_permres
    script:
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+-/,"")
       ourpheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
       base               =  bedfile.baseName
       phef               = "${base}_sh.phe"
       covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
       out                = "$our_pheno-$mperm"
       dir_gemma          =  "gemma"
       """
       shuf $phefi > $phef
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       gemma -bfile $base ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001 
       mkdir -p ${dir_gemma}/${this_pheno}
       mv output/* ${dir_gemma}/${this_pheno}/
       """
  }
  
  process doGemma{
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      set file(bedfile),file(bimfile),file(famfile) , file(rel_matrix),val(this_pheno), file(gemma_covariate), file(phefi) from gemma_data_perm2
    publishDir params.output_dir, overwrite:true, mode:'copy'
    output:
      file("${dir_gemma}/${out}.log.txt")
      set val(our_pheno), file("${dir_gemma}/${out}.assoc.txt") into gemma_manhatten_ch2
    script:
       our_pheno2         = this_pheno.replaceAll(/^[0-9]+-/,"")
       ourpheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
       our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
       base               =  bedfile.baseName
       covar_opt_gemma    =  (params.covariates) ?  " -c $gemma_covariate " : ""
       out                = "$our_pheno"
       dir_gemma          =  "gemma"
       """
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       gemma -bfile $base ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phefi -o $out -maf 0.0000001
       mkdir -p ${dir_gemma}
       mv output/* ${dir_gemma}/
       """
  }

// merge for each pheno process with gemma inital
gemma_merge=gemma_manhatten_ch2.join(gemma_permres.groupTuple())

process ComputePval{
   input :
     set val(pheno), file(gemmai), val(listgemperm) from gemma_merge
   publishDir params.output_dir, overwrite:true, mode:'copy'
   output :
     file("$out")
   script :
     tmplist=listgemperm.join("\n")
     filelist=pheno+"_listperm"
     out=gemmai+".newpval"
     """
     echo \"\"\"$tmplist\"\"\" >> $filelist
     perm_comp_pval.py --inp $gemmai --listgwas $filelist --out $out --head_pv "p_wald"
     """
}



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

nextflowversion =getres("nextflow -v")
if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"






