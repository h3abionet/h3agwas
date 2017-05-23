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
 *  2015-2016
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


/* Do permutation testing -- 0 for none, otherwise give number */
params.mperm = 1000

/* Adjust for multiple correcttion */
params.adjust = 0

supported_tests = ["chi2","fisher","model","cmh","linear","logistic"]


params.chi2     = 1
params.fisher   = 1
params.model    = 0
params.linear   = 0
params.logistic = 1
params.gemma    = 0
params.gemma_relopt = 1
params.gemma_lmmopt = 4


params.input_testing  = 'raw-GWA-data'

params.sexinfo_available = "false"


params.plink_process_memory = '750MB' // how much plink needs for this
params.other_process_memory = '750MB' // how much other processed need

max_plink_cores = params.max_plink_cores = 4

plink_mem_req = params.plink_process_memory
other_mem_req = params.other_process_memory

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




bed = Paths.get(params.input_dir,"${params.input_testing}.bed").toString()
bim = Paths.get(params.input_dir,"${params.input_testing}.bim").toString()
fam = Paths.get(params.input_dir,"${params.input_testing}.fam").toString()




pca_in_ch = Channel.create()
assoc_ch = Channel.create()
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .separate( pca_in_ch, assoc_ch ) { a -> [a,a] }





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

supported_tests = ["chi2","fisher","model","cmh","linear","logistic","gemma"]

requested_tests = supported_tests.findAll { entry -> params.get(entry)==1 }

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



pca_out_ch.subscribe { println "Drawn $it" }
/*comp_phase1_ch.subscribe { println "Done!!!" }*/
