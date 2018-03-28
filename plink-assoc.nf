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

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates","gemma_num_cores","gemma_mem_req","gemma","linear","logistic","chi2","fisher", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time"]



params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    println "Check $parm";
  }
}

def params_help = new LinkedHashMap(helps)


params.queue      = 'batch'
params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output_testing = "cleaned"
outfname = params.output_testing


/* Defines the path where any scripts to be executed can be found.
 */

println params.gemma_num_cores

/* Do permutation testing -- 0 for none, otherwise give number */
params.mperm = 1000

/* Adjust for multiple correcttion */
params.adjust = 0

supported_tests = ["chi2","fisher","model","cmh","linear","logistic"]


params.chi2     = 0
params.fisher   = 0
params.cmh     =  0
params.model   =  0
params.linear   = 0
params.logistic = 0
params.gemma = 1
params.gemma_mem_req = "6GB"
params.gemma_relopt = 1
params.gemma_lmmopt = 4


params.input_pat  = 'raw-GWA-data'

params.sexinfo_available = "false"


params.plink_mem_req = '750MB' // how much plink needs for this
params.other_process_memory = '750MB' // how much other processed need

max_plink_cores = params.max_plink_cores = 4

plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_memory

params.help = false


data_ch = Channel.fromPath(params.data)

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

//---- Modification of variables for pipeline -------------------------------//


def getConfig = {
  all_files = workflow.configFiles.unique()
  text = ""
  all_files.each { fname ->
      base = fname.baseName
      curr = "*-subsection{$base}@.@@.@*-footnotesize@.@*-begin{verbatim}"
      file(fname).eachLine { String line ->
	if (line.contains("secretKey")) { line = "secretKey='*******'" }
        if (line.contains("accessKey")) { line = "accessKey='*******'" }
        curr = curr + "@.@"+line 
      }
      curr = curr +"@.@*-end{verbatim}"
      text = text+curr
  }
  return text
}

/* Define the command to add for plink depending on whether sexinfo is
 * available or not.
 */
if ( params.sexinfo_available == "false" ) {
  sexinfo = "--allow-no-sex"
  println "Sexinfo not available, command --allow-no-sex\n"
} else {
  sexinfo = ""
  println "Sexinfo available command"
}



// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}




bed = Paths.get(params.input_dir,"${params.input_pat}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_pat}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_pat}.fam").toString()







gemma_assoc_ch = Channel.create()

pca_in_ch = Channel.create()
assoc_ch = Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .separate( pca_in_ch, assoc_ch, gemma_assoc_ch ) { a -> [a,a,a] }






process computePCA {
  cpus max_plink_cores
  memory plink_mem_req
  input:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from pca_in_ch

  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    set file("${outfname}.eigenval"), file("${outfname}.eigenvec")  \
         into pca_out_ch

  script:
  """
     plink --threads $max_plink_cores --bfile cleaned --pca --out ${outfname}
  """
}




num_assoc_cores = params.mperm == 0 ? 1 : 4

supported_tests = ["chi2","fisher","model","cmh","linear","logistic"]

requested_tests = supported_tests.findAll { entry -> params.get(entry) }


covariate = ""
gotcovar  = 0
pheno     = ""
if (params.covariates != "") {
    covariate = "--covar ${params.data} --covar-name ${params.covariates} "
    gotcovar = 1
}



 
if (params.data != "") {

  data_ch1 = Channel.create()
  data_ch2 = Channel.create()
  Channel.fromPath(params.data).separate(data_ch1,data_ch2) { a -> [a,a] } 
  
  process extractPheno {
    input:
     file(data) from data_ch1
    output:
     file(phenof) into pheno_ch
    script:
     phenof = "pheno.phe"
     """
     extractPheno.py $data ${params.pheno} $phenof
     """
  }


  pheno = "--pheno pheno.phe --all-pheno "




  process showPhenoDistrib {
    input:
    file(data) from data_ch2
    output:
      file ("B050*") into report_ch
    script:
      "phe_distrib.py --pheno ${params.pheno} $data B050 "
  }
}  else {
  report_ch = Channel.empty()
  pheno_ch  = Channel.from("dummy")
}



if (params.gemma == 1) {

  rel_ch = Channel.create()
  gem_ch = Channel.create()
  fam_ch = Channel.create()


  gemma_assoc_ch.separate (rel_ch, gem_ch, fam_ch) { a -> [a, a, a[2]] }

  process getGemmaRel {
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time params.big_time
    input:
       file plinks from rel_ch
    output:
       file("output/${base}.*XX.txt") into rel_mat_ch
    script:
       base = plinks[0].baseName
       """
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       gemma -bfile $base  -gk ${params.gemma_relopt} -o $base
       """
  }

   
  process transformCovariate {
    input:
      file(covariates) from data_ch 
      file(fam) from fam_ch
    output:
      set file(gemma_covariate), file(phef) into gemma_data_ch
      stdout into pheno_cols_ch
    script:
      base = fam.baseName
      gemma_covariate = "${base}.gemma_cov"
      phef = "${base}_n.phe"
      """
      gemma_covariate.py --data  $covariates --inp_fam  $fam --cov_list ${params.covariates} \
                          --pheno ${params.pheno} --cov_out $gemma_covariate --phe_out ${phef}
      """
  }

  ind_pheno_cols_ch = Channel.create()
  check = Channel.create()
  pheno_cols_ch.flatMap { list_str -> list_str.split() }.tap ( check) .set { ind_pheno_cols_ch }

  check.subscribe { println "Found valuye phe $it" }

  process doGemma {
    cpus 2
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(plinks) from  gem_ch
      file(matrix) from  (rel_mat_ch)
      set file (covariate), file (phef) from gemma_data_ch
    each this_pheno from ind_pheno_cols_ch
    publishDir params.output_dir
    output:
      file("output/${base}-${this_pheno}.log.txt")
      set val(base), val(this_pheno), file("output/${base}-${this_pheno}.assoc.txt") into gemma_manhatten_ch
    script:
      base = plinks[0].baseName
      """
      export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
      this_pheno_col=`echo ${this_pheno} | sed "s/-.*//"`
      gemma -bfile $base  -c $covariate  -k $matrix -lmm 1  -n \${this_pheno_col} -p $phef -o $base-$this_pheno
      """
  }



  process showGemmaManhatten { 
    publishDir params.output_dir
    input:
      set val(base), val(this_pheno), file(assoc) from gemma_manhatten_ch
    output:
      set file("${out}*")  into report_gemma_ch
    script:
      out = "C049$this_pheno"
      """
      gemma_man.py  $assoc $this_pheno ${out}
      """
  }

  report_ch = report_ch.flatten().mix(report_gemma_ch.flatten())
    
} 

    


if (params.chi2+params.fisher+params.logistic+params.linear > 0) {

   process computeTest {
      cpus num_assoc_cores
      input:
       set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from assoc_ch    
       file (phenof) from pheno_ch
      each test from requested_tests
      publishDir params.output_dir, overwrite:true, mode:'copy'
      output:
      set val(test), file("${outfname}.*") into out_ch
      script:
       base = "cleaned"
       perm = (params.mperm == 0 ? "" : "mperm=${params.mperm}")
       adjust = (params.adjust ? "--adjust" : "")
       template "${test}.sh"
   }



  process drawPlinkResults { 
    input:
      set val(test), file(results) from out_ch
    output:
      set file("${base}*man*png"), file ("${base}*qq*png"), file("C050*tex") into report_plink
    publishDir params.output_dir
    script:
      base="cleaned"
      """
      plinkDraw.py  $base $test "${params.pheno}" $gotcovar png
      """
  }

  report_ch = report_ch.mix(report_plink.flatten())
  
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


process doReport {
  input:
    file(reports) from report_ch.toList()
  publishDir params.output_dir
  output:
    file("${out}.pdf")
  script:
    out = params.output+"-report"
    config = getConfig()
    images = workflow.container
    texf   = "${out}.tex"
    template "make_assoc_report.py"
}




pca_out_ch.subscribe { println "Drawn $it" }
