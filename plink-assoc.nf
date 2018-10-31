#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *
 *      Scott Hazelhurst
 *      Jean-Tristan Brandenburg
 *      Shaun Aron
 *   	Rob Clucas
 *      Eugene de Beste
 *      Lerato Magosi
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2018
 *
 *(C) University of the Witwatersrand, Johannesburg, 2016-2018 on behalf of the H3ABioNet Consortium
 *This is licensed under the MIT Licence. See the "LICENSE" file for details
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths


def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates","gemma_num_cores","gemma_mem_req","gemma","linear","logistic","chi2","fisher", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "super_pi_hat", "f_lo_male", "f_hi_female", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "AMI", "instanceType", "instance-type", "bootStorageSize", "boot-storage-size", "maxInstances", "max-instances", "other_mem_req", "sharedStorageMount", "shared-storage-mount", "max_plink_cores", "pheno","big_time","thin","adjust","mperm"]



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



/* Defines the path where any scripts to be executed can be found.
 */



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
params.gemma = 0
params.gemma_mem_req = "6GB"
params.gemma_relopt = 1
params.gemma_lmmopt = 4


params.input_pat  = 'raw-GWA-data'

params.sexinfo_available = "false"


params.plink_mem_req = '750MB' // how much plink needs for this
params.other_mem_req = '750MB' // how much other processed need

max_plink_cores = params.max_plink_cores = 4

plink_mem_req = params.plink_mem_req
other_mem_req = params.other_mem_req

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

pca_in_ch = Channel.create()
prune_in_ch = Channel.create()
assoc_ch  = Channel.create()
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
if (params.chi2) println "Doing chi2 testing"
if (params.linear) println "Doing linear regression testing"
if (params.logistic) println "Doing logistic regression testing"
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
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  ( prune_in_ch, pca_in_ch, assoc_ch, gemma_assoc_ch )
    script:
       base = bed.baseName
       out  = base+"_t"
       "plink --bfile $base $thin $chrom --make-bed --out $out"
  }

  println "\nData has been thinned or only some chromosomes used  (is the run for test purposes only?)\n"
   


} else {
  raw_src_ch.separate( prune_in_ch, pca_in_ch, assoc_ch, gemma_assoc_ch ) { a -> [a,a,a,a] }
}



process pruneData {
  cpus 1
  memory plink_mem_req
  input:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from prune_in_ch
  output:
    file("check.prune.in") into list_prune_ch
  script:
      base = "cleaned"
     """
     plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
     """
}


process computePCA {
  cpus max_plink_cores
  memory plink_mem_req
  input:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from pca_in_ch
    file("check.prune.in") from list_prune_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    set file("${outfname}.eigenval"), file("${outfname}.eigenvec")  \
         into pca_out_ch
  script:
      base = "cleaned"
      prune= "${base}-prune"
     """
     plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
     plink --threads $max_plink_cores --bfile $prune --pca --out ${outfname}
     """
}

process drawPCA {
    input:
      set file(eigvals), file(eigvecs) from pca_out_ch
    output:
      file (output) into report_pca_ch
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





num_assoc_cores = params.mperm == 0 ? 1 : Math.min(10,params.max_plink_cores)

supported_tests = ["chi2","fisher","model","cmh","linear","logistic"]

requested_tests = supported_tests.findAll { entry -> params.get(entry) }


covariate = ""
gotcovar  = 0
pheno     = ""



 
if (params.data != "") {

   checker(file(params.data))

   if (params.covariates != "") {
      gotcovar = 1
  }

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
     all_phenos = params.covariates.length()>0 ? params.pheno+","+params.covariates : params.pheno
     """
     extractPheno.py $data ${all_phenos} $phenof
     """
  }


  pheno_label_ch = Channel.from(params.pheno.split(","))

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
  pheno_label = ""
  pheno_label_ch = Channel.from("")
}





if (params.gemma == 1) {

  rel_ch = Channel.create()
  gem_ch = Channel.create()
  fam_ch = Channel.create()


  gemma_assoc_ch.separate (rel_ch, gem_ch, fam_ch) { a -> [a, a, a[2]] }

  process getGemmaRel {
    label 'gemma'
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

   
  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""
  
  process  getGemmaPhenosCovar {
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
      gemma_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
                          --pheno ${params.pheno} --cov_out $gemma_covariate --phe_out ${phef}
      """
  }

  ind_pheno_cols_ch = Channel.create()
  check = Channel.create()
  pheno_cols_ch.flatMap { list_str -> list_str.split() }.tap ( check) .set { ind_pheno_cols_ch }

  check.subscribe { println "Found phenotype request $it" }

  process doGemma {
    label 'gemma'
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(plinks) from  gem_ch
      file(matrix) from  (rel_mat_ch)
      set file (covariate), file (phef) from gemma_data_ch
    each this_pheno from ind_pheno_cols_ch
    publishDir params.output_dir
    output:
      file("gemma/${out}.log.txt")
      set val(base), val(our_pheno), file("gemma/${out}.assoc.txt") into gemma_manhatten_ch
    script:
      base = plinks[0].baseName
      covar_opt =  (params.covariates) ?  " -c $covariate" : ""
      our_pheno = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
      out = "$base-$our_pheno"
      """
      export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
      this_pheno_col=`echo ${this_pheno} | sed 's/-.*//' `
      gemma -bfile $base ${covar_opt}  -k $matrix -lmm 1  -n \${this_pheno_col} -p $phef -o $out
      mv output gemma
      """
  }



  process showGemmaManhatten { 
    publishDir params.output_dir
    input:
      set val(base), val(this_pheno), file(assoc) from gemma_manhatten_ch
    output:
      file("${out}*")  into report_gemma_ch
    script:
      our_pheno = this_pheno.replaceAll("_","-")
      out = "C049$this_pheno"
      """
      gemma_man.py  $assoc $this_pheno ${out}cd 
      """
  }

  report_ch = report_ch.flatten().mix(report_gemma_ch.flatten())
    
} 


    


if (params.chi2+params.fisher+params.logistic+params.linear > 0) {

   process computeTest {
      cpus num_assoc_cores
      time params.big_time
      input:
       set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from assoc_ch    
       file (phenof) from pheno_ch
      each test_choice from requested_tests
      each pheno_name from pheno_label_ch
      publishDir "${params.output_dir}/${test}", overwrite:true, mode:'copy'
      output:
        set val(test), val(pheno_name), file("${outfname}.*") into out_ch
      script:
       base = "cleaned"
       pheno_name = pheno_name.replaceFirst("/.*","")
       perm = (params.mperm == 0 ? "" : "mperm=${params.mperm}")
       adjust = (params.adjust ? "--adjust" : "")
       outfname = "${pheno_name}"
       test = test_choice == "chi2" ? "assoc" : test_choice
       if (params.data == "") {
           pheno_cmd = ""
           out = base
       } else {
           pheno_cmd = "--pheno $phenof --pheno-name $pheno_name "
           if (params.covariates) covariate = "--covar ${phenof} --covar-name ${params.covariates} "
           out = pheno
       }
       template "test.sh"
   }


  log_out_ch = Channel.create()
 
  log_out_ch.subscribe { println "Completed plink test ${it[0]}" }
 
  process drawPlinkResults { 
    input:
    set val(test), val(pheno_name), file(results) from out_ch.tap(log_out_ch)
    output:
      set file("${base}*man*pdf"), file ("${base}*qq*pdf"), file("C050*tex") into report_plink
    publishDir params.output_dir
    script:
      base="cleaned"
      """
      plinkDraw.py  C050 $base $test ${pheno_name} $gotcovar pdf
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

report_ch = report_ch.mix(report_pca_ch)

process doReport {
  label 'latex'
  input:
    file(reports) from report_ch.toList()
  publishDir params.output_dir
  output:
    file("${out}.pdf") into final_report_ch
  script:
    out = params.output+"-report"
    these_phenos     = params.pheno
    these_covariates = params.covariates
    config = getConfig()
    images = workflow.container
    texf   = "${out}.tex"
    template "make_assoc_report.py"
}


final_report_ch.subscribe { 
     b=it.baseName; println "The output report is called ${params.output_dir}/${b}.pdf"
     params.pheno.split(",").each { p ->
       if (covs.contains(p)) {
         println("\n\nThe phenotype <$p> is also given as a covariate -- this seems like a very bad idea")
       }
     }
}


