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


/* Defines the path where any scripts to be executed can be found.
 */

params.scripts   = "${params.work_dir}/scripts"


       /* What association tests should be done
	*/




/* Do permutation testing -- 0 for none, otherwise give number */
params.mperm = 1000

/* Adjust for multiple correcttion */
params.adjust = true

supported_tests = ["chi2","fisher","model","cmh","linear","logistic"]
params.assoc = ["chi2","fisher","model","cmh","linear","logistic"]



/* Defines the names of the plink binary files in the plink directory
 * (.fam, .bed, .bed).
 *
 * NOTE: This must be without the extension (so if A.fam, A.bed, ...
 *       then use 'A').
 */
params.data_name  = 'raw-GWA-data'

/* When computing IBD do we want to exclude high-lD regions from computation */
/* empty string if not */

params.high_ld_regions_fname = ""

/* Defines if sexinfo is available or not, options are:
 *  - "true"  : sexinfo is available
 *  - "false" : sexinfo is not avalable
 */
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

// From the input base file, we get the bed, bim and fam files -- absolute path and add suffix

bim = Paths.get(params.input_dir,"${params.plink_fname}.bim").toString()
//println(bim)
//println(input_dir)

// Prepends scripts directory path to the argument given
def path = {
   fn ->
     Paths.get(params.scripts,fn).toString()
}

// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}


//------------

/* Deal with scripts -- we check if the scripts exist and if they do */




all_scripts.each  { checker(file(path(it))) }

script_ch = [ 'dummy' : Channel.empty()]

all_scripts.each { 
   name-> 
     script_ch.put(name, Channel.fromPath(path(name)))
}




// Creating two channels with the file names and at the same time
// checking file existence





bed = Paths.get(params.input_dir,"${params.data_name}.bed").toString()
bim = Paths.get(params.input_dir,"${params.data_name}.bim").toString()
fam = Paths.get(params.input_dir,"${params.data_name}.fam").toString()




bim_ch = Channel.fromPath(bim).map checker
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
    set file('cleaned.eigenval'), file('cleaned.eigenvec')  into pca_out_ch

  script:
  """
     plink --threads $max_plink_cores --bfile cleaned --pca --out cleaned
  """
}



num_assoc_cores = params.mperm == 0 ? 1 : max_plink_cores

process computeTest {
   echo true
   cpus num_assoc_cores
   input:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from assoc_ch    
   each test from params.assoc
   publishDir params.output_dir, overwrite:true, mode:'copy'
   output:
      set file("cleaned.*") into out_ch
   script:
    base = "cleaned"
    perm = (params.mperm == 0 ? "" : "mperm=${params.mperm}")
    adjust = (params.adjust ? "--adjust" : "")
    if (test == "chi2")
      template "chi2.sh"
    else if (test == "fisher") 
      template "fisher.sh"
}



pictures_ch.subscribe { println "Drawn $it" }
/*comp_phase1_ch.subscribe { println "Done!!!" }*/
