#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *
 *      Shaun Aron
 *   	Rob Clucas
 *      Eugene de Beste
 *      Scott Hazelhurst
 *      Anmol Kiran
 *      Lerato Magosi
 *      Abayomi Mosaku
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2017
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths


def helps = [ 'help' : 'help' ]


def params_help = new LinkedHashMap(helps)


params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"
params.output_dir = "${params.work_dir}/output"
params.output_testing = "cleaned"
outfname = params.output_testing


/* Defines the path where any scripts to be executed can be found.
 */

println params.gemma_plink_cores

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


params.plink_process_memory = '750MB' // how much plink needs for this
params.other_process_memory = '750MB' // how much other processed need

max_plink_cores = params.max_plink_cores = 4

plink_mem_req = params.plink_process_memory
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




num_assoc_cores = params.mperm == 0 ? 1 : max_plink_cores

supported_tests = ["chi2","fisher","model","cmh","linear","logistic"]

requested_tests = supported_tests.findAll { entry -> params.get(entry) }


covariate = ""
pheno     = ""
if (params.covariate != "") {
    covariate = "--covar ${params.data} --covar-name ${params.covariate} "
}
 
if (params.data != "") {
    pheno = "--pheno ${params.data} --pheno-name ${params.pheno}"
}


process computeTest {
   echo true
   cpus num_assoc_cores
   input:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from assoc_ch    
   each test from requested_tests
   publishDir params.output_dir, overwrite:true, mode:'copy'
   output:
      set file("${outfname}.*") into out_ch
   script:
    base = "cleaned"
    perm = (params.mperm == 0 ? "" : "mperm=${params.mperm}")
    adjust = (params.adjust ? "--adjust" : "")
    template "${test}.sh"
}


if (params.gemma == 1) {

  rel_ch = Channel.create()
  gem_ch = Channel.create()
  fam_ch = Channel.create()


  gemma_assoc_ch.separate (rel_ch, gem_ch, fam_ch) { a -> [a, a, a[2]] }

  process getGemmaRel {
    cpus params.gemma_plink_cores
    memory params.plink_mem_req
    input:
       file plinks from rel_ch
    output:
       file("output/$base}.sXX.txt") into rel_mat_ch
    script:
       base = plinks[0].baseName
       """
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
      export OPENBLAS_NUM_THREADS=${params.gemma_plink_cores}
      gemma_covariate.py --data  $covariates --inp_fam  $fam --cov_list ${params.covariates} \
                          --pheno ${params.pheno} --cov_out $gemma_covariate --phe_out ${phef}
      """
  }

  process doGemma {
    cpus params.gemma_plink_cores
    memory params.gemma_mem_req
    input:
      file(plinks) from  gem_ch
      file(matrix) from  rel_mat_ch
      set file (covariate), file (new_fam) from gemma_data_ch
      val(pheno_cols) from pheno_cols_ch
    publishDir params.output_dir
    output:
      file("${base}.log.txt")
      set val(base), file("${base}.assoc.txt") into gemma_manhatten_ch
    script:
      base = plinks[0]
      """
      export OPENBLAS_NUM_THREADS=${params.gemma_plink_cores}
      gemma -bfile $base  -c $covariate  -k $matrix -lmm 1  -n ${pheno_cols}  -o $base
      """
  }

  process showGemmaManhatten { 
    publishDir params.output_dir
    input:
      set val(base), file(assoc) from gemma_manhatten_ch
    output:
      file(out)  
    script:
      out = "${base}.png"
      """
      gemma_man.py  $assoc $png
      """
  }
    
} 


if (params.chi2+params.fisher+params.logistic+params.linear > 0) {

  process drawPlinkResults { 
    input:
      file(results) from out_ch
    output:
      set file(man), file (qq), file(tex) into report_ch
    publishDir params.output_dir
    script:
      base=results[0].baseName
      man ="${base}-man.png"
      qq  ="${base}-qq.png"
      tex ="C050.tex"
      """
      plinkDraw.py  $man $qq $tex
      """
  }

}

pca_out_ch.subscribe { println "Drawn $it" }



