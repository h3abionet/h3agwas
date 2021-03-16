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

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths




def helps = [ 'help' : 'help' ]

allowed_params = ["input_dir","input_pat","output","output_dir","data","plink_mem_req","covariates","gemma_num_cores","gemma_mem_req","gemma","linear","logistic","assoc","fisher", "work_dir", "scripts", "max_forks", "high_ld_regions_fname", "sexinfo_available", "cut_het_high", "cut_het_low", "cut_diff_miss", "cut_maf", "cut_mind", "cut_geno", "cut_hwe", "pi_hat", "case_control", "case_control_col", "phenotype", "pheno_col", "batch", "batch_col", "samplesize", "strandreport", "manifest", "idpat", "accessKey", "access-key", "secretKey", "secret-key", "region", "other_mem_req", "max_plink_cores", "pheno","big_time","thin", "gemma_mat_rel","print_pca", "file_rs_buildrelat","genetic_map_file", "rs_list","adjust","bootStorageSize","shared-storage-mount","mperm","sharedStorageMount","max-instances","maxInstances","boot-storage-size","sharedStorageMound","instance-type","instanceType","AMI", "gemma_multi", "sample_snps_rel"]


/*bolt_use_missing_cov --covarUseMissingIndic : “missing indicator method” (via the --covarUseMissingIndic option), which adds indicator variables demarcating missing status as additional covariates. */
ParamBolt=["bolt_ld_scores_col", "bolt_ld_score_file","boltlmm", "bolt_covariates_type",  "bolt_use_missing_cov", "bolt_num_cores", "bolt_mem_req", "exclude_snps", "bolt_impute2filelist", "bolt_impute2fidiid", "bolt_otheropt","bolt_bin"]
allowed_params+=ParamBolt
ParamFast=["fastlmm","fastlmm_num_cores", "fastlmm_mem_req", "fastlmm_multi", "fastlmmc_bin", "covariates_type"]
allowed_params+=ParamFast
/*Gxe : */
GxE_params=['gemma_gxe', "plink_gxe", "gxe"]
allowed_params+=GxE_params

FastGWA_params=["fastgwa_memory", "fastgwa_cpus", 'fastgwa']
allowed_params+=FastGWA_params

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
params.mperm = 00

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
params.gemma_multi=0
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

params.fastlmmc_bin ="fastlmmc"
params.bolt_bin ="bolt"
params.gemma_bin ="gemma"

/*gxe param : contains column of gxe*/
params.gemma_gxe=0
params.plink_gxe=0
params.max_plink_cores = 4
params.rs_list=""
params.gxe=""

/**/
params.fastgwa=0
params.fastgwa_memory="10G"
params.fastgwa_cpus=5
params.grm_nbpart=100
params.gcta64_bin = "gcta64"
params.fastgwa_type="--fastGWA-mlm-exact"
params.grm_cutoff =  0.05
params.covariates_type=""
params.gcta_grmfile=""
params.sample_snps_rel=0


params.input_pat  = 'raw-GWA-data'

params.sexinfo_available = "false"


params.plink_mem_req = '6GB' // how much plink needs for this
params.other_process_memory = '10G' // how much other processed need


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
assoc_ch_gxe_freq=Channel.create()
grlm_assoc_ch = Channel.create()
ch_bolt_snpchoice=Channel.create()
fastgwa_assoc_ch = Channel.create()
raw_src_ch= Channel.create()
ch_select_rs_format=Channel.create()

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
if (params.gemma) println "Doing gemma testing"
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
      set file("${out}.bed"), file("${out}.bim"), file("${out}.fam") into  ( pca_in_ch, assoc_ch, gemma_assoc_ch, boltlmm_assoc_ch,fastlmm_assoc_ch, rel_ch_fastlmm,assoc_ch_gxe_freq, assoc_ch_gxe, grlm_assoc_ch, fastgwa_assoc_ch, ch_bolt_snpchoice, ch_select_rs_format)
    script:
       base = bed.baseName
       out  = base+"_t"
       "plink --keep-allele-order --bfile $base $thin $chrom --make-bed --out $out"
  }

  println "\nData has been thinned or only some chromosomes used  (is the run for test purposes only?)\n"
   


} else {
    /*JT : append boltlmm_assoc_ch and a]*/
    raw_src_ch.separate( pca_in_ch, assoc_ch, assoc_ch_gxe,assoc_ch_gxe_freq,gemma_assoc_ch, boltlmm_assoc_ch, fastlmm_assoc_ch,rel_ch_fastlmm, grlm_assoc_ch,fastgwa_assoc_ch,ch_bolt_snpchoice,ch_select_rs_format) { a -> [a,a,a,a,a,a,a,a,a,a,a, a,a] }
}




if(params.print_pca!=0){
   process computePCA {
     cpus max_plink_cores
     memory plink_mem_req
     time   params.big_time
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
	plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
	plink --keep-allele-order --bfile ${base} --extract check.prune.in --make-bed --out $prune
	plink --threads $max_plink_cores --bfile $prune --pca --out ${outfname}
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





num_assoc_cores = params.mperm == 0 ? 1 : Math.min(10,params.max_plink_cores)

supported_tests = ["assoc","fisher","model","cmh","linear","logistic"]

requested_tests = supported_tests.findAll { entry -> params.get(entry) }


covariate = ""
gotcovar  = 0
pheno     = ""

/*Case where we sample builrelatdness*/
balise_filers_rel=1
if(params.boltlmm+params.gemma+params.fastlmm+params.fastgwa>0){
 if(params.file_rs_buildrelat=="" && params.sample_snps_rel==1){
  process select_rs_format{
     cpus max_plink_cores
     memory plink_mem_req
     time   params.big_time
     input :
       set file(bed),file(bim), file(fam) from ch_select_rs_format
    output:
       file("${prune}.prune.in") into  filers_matrel_mat_fast, filers_matrel_mat_GWA, filers_matrel_mat_gem, filers_matrel_bolt, filers_count_line
     script:
        base = bed.baseName
        prune= "${base}-prune"
        """
        plink --bfile ${base} --indep-pairwise 100 20 0.1 --out $prune
        """
   }
   //BoltNbMaxSnps=filers_count_line.countLines()
 }else{
/* n*/
 if(params.file_rs_buildrelat==""){
   balise_filers_rel=0
   filers_matrel_mat_fast=file('NO_FILE')
   filers_matrel_mat_GWA=file('NO_FILE')
   filers_matrel_mat_gem=file('NO_FILE')
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
  }
 }
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
    // not sure difference between container and label
    input:
    file(data) from data_ch
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
       label 'py2fast'
       cpus params.fastlmm_num_cores
       time   params.big_time
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
	 our_pheno = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/-$/,"")
	 covar_opt_fast =  (params.covariates) ?  " -covar $covariate" : ""
	 newbase=base+"-"+chro
	 out = "$base-$our_pheno"+"-"+chro+".stat"
	 """
	 this_pheno_col=`echo ${this_pheno} | sed 's/-.*//' `
	 plink --keep-allele-order --bfile $base --chr $chro --make-bed --out $newbase
	 $fastlmmc -REML -simType RRM -verboseOut -sim $rel -bfile $newbase -pheno ${phef} -simLearnType Full -out $out -maxThreads $params.fastlmm_num_cores \
	          $covar_opt_fast  -mpheno \${this_pheno_col} -bfileSim $base
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
       label 'py2fast'
       cpus params.fastlmm_num_cores
       time   params.big_time
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
	 this_pheno_col=`echo ${this_pheno} | sed 's/-.*//' `
	 $fastlmmc -REML -simType RRM -verboseOut -bfile $base -pheno ${phef} -simLearnType Full -out $out -maxThreads $params.fastlmm_num_cores \
	           $covar_opt_fast  -mpheno \${this_pheno_col} -bfileSim $base
	 """
     }
  }

  // this part is plotting done for any fastlmm mode
  //, overwrite:true, mode:'copy'
  process showFastLmmManhatten {
    memory params.other_process_memory
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
  process doCountNbSnp{
    time   params.big_time
    input :
       file(bim) from bim_ch_bolt
    output :
       stdout into nbsnp_ch_bolt
    script :
      """
      wc -l $bim|awk '{print \$1}'
      """
  }
  /*    nb_snp= CountLinesFile(base+".bim") */
  if(params.exclude_snps)rs_ch_exclude_bolt=Channel.fromPath(params.exclude_snps, checkIfExists:true)
  else rs_ch_exclude_bolt=file('NO_FILE')

  if(params.bolt_impute2filelist!=""){
  Impute2FileList=Channel.fromPath(params.bolt_impute2filelist, checkIfExists:true)
  Impute2FID = Channel.fromPath(params.bolt_impute2fidiid, checkIfExists:true)
  }else{
  Impute2FileList=file('NO_FILE1')
  Impute2FID = file('NO_FILE2')
  }
  if(params.bolt_ld_score_file!=""){
     Bolt_ld_score= Channel.fromPath(params.bolt_ld_score_file, checkIfExists:true)
  }else{
     Bolt_ld_score = file('NO_FILE3')
  }
//genetic_map_file
  if(params.genetic_map_file!=""){
     Bolt_genetic_map= Channel.fromPath(params.genetic_map_file, checkIfExists:true)
  }else{
     Bolt_genetic_map = file('NO_FILE4')
  }

  process doBoltmm{
    label 'py2fast'
    cpus params.bolt_num_cores
    memory params.bolt_mem_req
    time   params.big_time
    input:
      set file(plinksbed), file(plinksbim), file(plinksfam) from plink_ch_bolt
      val nb_snp from nbsnp_ch_bolt
      file(phef) from newdata_ch_bolt
      file(rs_exclude) from rs_ch_exclude_bolt
      file(SnpChoiceMod) from filers_matrel_bolt
      file(imp2_filelist) from Impute2FileList 
      file(imp2_fid) from Impute2FID
      file(bolt_ld_score) from Bolt_ld_score
      file(bolt_genetic_map) from Bolt_genetic_map
    publishDir "${params.output_dir}/boltlmm", overwrite:true, mode:'copy'
    each this_pheno from ind_pheno_cols_ch_bolt
    output:
      file(outbolt)
      set val(base), val(our_pheno), file("$outf") into bolt_manhatten_ch
    script:
      base = plinksbed.baseName
      our_pheno2         = this_pheno.replaceAll(/^[0-9]+@@@/,"")
      our_pheno3         = our_pheno2.replaceAll(/\/np.\w+/,"")
      our_pheno          = this_pheno.replaceAll(/_|\/np.\w+/,"-").replaceAll(/[0-9]+@@@/,"")
      outimp  = (params.bolt_impute2filelist!="") ? "$base-${our_pheno2}.imp.stat" : "$base-${our_pheno2}.stat"
      outbolt     = "$base-${our_pheno2}.stat" 
      outf    = (params.bolt_impute2filelist!="") ? outimp : outbolt
      outReml = "$base-$our_pheno2"+".reml"
      covar_file_bolt =  (params.covariates) ?  " --covarFile ${phef} " : ""
      model_snp  = "--modelSnps=$SnpChoiceMod "
      ld_score_cmd = (params.bolt_ld_score_file!="") ? "--LDscoresFile=$bolt_ld_score" :" --LDscoresUseChip "
      ld_score_cmd = (params.bolt_ld_score_file!="" & params.bolt_ld_scores_col!="") ? "$ld_score_cmd --LDscoresCol=${params.bolt_ld_scores_col}" :" $ld_score_cmd "
      exclude_snp = (params.exclude_snps!="") ? " --exclude $rs_exclude " : ""
      boltimpute = (params.bolt_impute2filelist!="") ? " --impute2FileList $imp2_filelist --impute2FidIidFile $imp2_fid --statsFileImpute2Snps $outimp  " : ""
      geneticmap = (params.genetic_map_file!="") ?  " --geneticMapFile=$bolt_genetic_map " : ""
      """
      BoltNbMaxSnps=`cat  ${SnpChoiceMod}|wc -l`
      bolt.py ${params.bolt_bin} $type_lmm --bfile=$base  --phenoFile=${phef} --phenoCol=${our_pheno3} \
     --numThreads=$params.bolt_num_cores $cov_bolt $covar_file_bolt --statsFile=$outbolt \
    $ld_score_cmd  $missing_cov --lmmForceNonInf  $model_snp $exclude_snp $boltimpute $geneticmap ${params.bolt_otheropt} \
      --maxModelSnps=\$BoltNbMaxSnps
      #bolt.py bolt  --reml  --bfile=$base  --phenoFile=${phef} --phenoCol=${our_pheno3} --numThreads=$params.bolt_num_cores $cov_bolt $covar_file_bolt $missing_cov $model_snp $geneticmap $exclude_snp |\
      #       grep -B 1 -E "^[ ]+h2" 1> $outReml 
      """
  }

  process showBoltmmManhatten {
   memory params.other_process_memory
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


}/*JT End of boltlmm*/ else {
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

if(params.gemma_multi==1){
     process getListeChroGem{
        input :
          file(BimFile) from bim_ch_fast_gem
        output :
          stdout into (chrolist,chrolist2, chrolisti_gxe)
        script:
         """
         cat $BimFile|awk '{print \$1}'|uniq|sort|uniq
        """
     }
  check2 = Channel.create()
  list_chro_gemma=chrolist.flatMap { list_str -> list_str.split() }.tap ( check2)
  check2 = Channel.create()
  list_chro_gemma_gxe=chrolisti_gxe.flatMap { list_str -> list_str.split() }.tap ( check2)



}
if (params.gemma == 1){
  ind_pheno_cols_ch = newNamePheno(params.pheno)

  if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""
   if(params.rs_list==""){
        rsfile=file('NO_FILE5')
     }else{
        rsfile=file(params.rs_list)
   }


 if(params.gemma_multi==1){

  check2 = Channel.create()

  process doGemmaChro{
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(covariates) from data_ch
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
       list_ind_nomissing.py --data $covariates --inp_fam $inp_fam $covariate_option --pheno $our_pheno3 --dataout $data_nomissing \
                             --lindout $list_ind_nomissing
       gemma_relselind.py  --rel $rel --inp_fam $inp_fam --relout $rel_matrix --lind $list_ind_nomissing
       plink --keep-allele-order --bfile $base --keep $list_ind_nomissing --make-bed --out $newbase  ${rs_plk_gem} --chr $chro
       all_covariate.py --data  $data_nomissing --inp_fam  ${newbase}.fam $covariate_option --cov_out $gemma_covariate \
                          --pheno $our_pheno2 --phe_out ${phef} --form_out 1
       export OPENBLAS_NUM_THREADS=${params.gemma_num_cores}
       ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001
       mv output ${dir_gemma}
       rm $rel_matrix
       rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
       """
}
     gemma_manhatten_ch_chro_merge=gemma_manhatten_ch_chro.groupTuple()
     process doMergeGemma{
          input :
            set (val(this_pheno),list_file, base_list) from  gemma_manhatten_ch_chro_merge
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

  }else{
  process doGemma{
    cpus params.gemma_num_cores
    memory params.gemma_mem_req
    time   params.big_time
    input:
      file(covariates) from data_ch
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
       ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001 
       mv output ${dir_gemma}
       rm $rel_matrix
       rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
       """
  }
  }


  process showGemmaManhatten {
    memory params.other_process_memory
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
        rsfile=file('NO_FILERS')
     }else{
        rsfile=file(params.rs_list)
   }

 if(params.gemma_multi==1){
  ind_pheno_cols_ch_gxe_multi = newNamePheno(params.pheno)
process doGemmaGxEChro{
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
    publishDir params.output_dir, overwrite:true, mode:'copy'
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
       ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001 $gxe_opt_gemma
       mv output ${dir_gemma}
       rm ${newbase}.bed ${newbase}.bim ${newbase}.fam
       """
  }


     gemma_manhatten_ch_chro_gxe_merge=gemma_manhatten_ch_chro_gxe.groupTuple()
     process doMergeGemmaGxE{
          input :
            set (val(this_pheno),list_file, base_list) from  gemma_manhatten_ch_chro_gxe_merge
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
  }else{

   ind_pheno_cols_ch = newNamePheno(params.pheno)

   process doGemmaGxE{
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
       ${params.gemma_bin} -bfile $newbase ${covar_opt_gemma}  -k $rel_matrix -lmm 1  -n 1 -p $phef -o $out -maf 0.0000001 $gxe_opt_gemma
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
       added_freq_gxe.py --bfile $base --file_gxe $gemmares --pheno_file $data_ch --pheno ${our_pheno3} --pheno_gxe ${params.gxe} --out $gemmaresfreq  --plk_cores ${params.gemma_num_cores}
       """
  } 

  process showGemmaManhattenGxE { 
    memory params.other_process_memory
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
       template "test.sh"
   }


 // log_out_ch = Channel.create()
 
  //log_out_ch.subscribe { println "Completed plink test ${it[0]}" }
 
  process drawPlinkResults { 
    memory params.other_process_memory
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
   if(params.rs_list==""){
        rsfile_plkgxe=file('NO_FILERS')
     }else{
        rsfile_plkgxe=file(params.rs_list, checkIfExists=true)
   }


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

  process showPlinkManhattenGxE {
    memory params.other_process_memory
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

// Partition the GRM into 100 parts
// ${params.gcta64_bin} --bfile test --make-grm-part 100 1 --thread-num 5 --out test
// ${params.gcta64_bin} --bfile test --make-grm-part 100 2 --thread-num 5 --out test
// ...
// ${params.gcta64_bin} --bfile test --make-grm-part 100 100 --thread-num 5 --out test
// # Merge all the parts together (Linux, Mac)
// cat test.part_3_*.grm.id > test.grm.id
// cat test.part_3_*.grm.bin > test.grm.bin
// cat test.part_3_*.grm.N.bin > test.grm.N.bin
// # Make a sparse GRM from the merged full-dense GRM
// ${params.gcta64_bin} --grm test_grm --make-bK-sparse 0.05 --out test_sp_grm
//
if(params.fastgwa==1){

if(params.gcta_grmfile==""){
 process FastGWADoGRM{
    memory params.fastgwa_memory
    cpus params.fastgwa_cpus
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
     ${params.gcta64_bin} --bfile $base --make-grm-part  ${params.grm_nbpart} $mpart --thread-num ${params.fastgwa_cpus} --out mgrm $rs_list
     """
 }


  idgrm_c=idgrm.collect()
  bingrm_c=bingrm.collect()
  nbingrm_c=nbingrm.collect()
  process MergFastGWADoGRM{
    memory params.fastgwa_memory
    cpus params.fastgwa_cpus
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
      ${params.gcta64_bin} --grm test_grm --make-bK-sparse ${params.grm_cutoff} --out $head --thread-num ${params.fastgwa_cpus}
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

balqualcov=params.covariates_type!="" & params.covariates_type.split(',').contains('1') 
balquantcov=params.covariates_type!="" & params.covariates_type.split(',').contains('0') 
process FastGWARun{
    memory params.fastgwa_memory
    cpus params.fastgwa_cpus
    input :
       set val(head),file(alldigrm), file(allbingrm) from grm_all
       set file(bed),file(bim),file(fam) from fastgwa_assoc_ch 
       file(covariates) from data_ch_fastgwa
    publishDir "${params.output_dir}/fastgwa/", overwrite:true, mode:'copy'
    each this_pheno from pheno_spl_gcta
    output :
       set val(base), val(this_pheno), file("${out}.fastGWA") into fastgwa_manhatten_ch
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
     """
     all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option --pheno ${this_pheno} --phe_out ${phef}  --cov_out $covfilequant --form_out 4  $covqual_cov
     ${params.gcta64_bin} --bfile $base ${params.fastgwa_type}  --pheno $phef  $covquant_fastgwa --threads ${params.fastgwa_cpus} --out $out --grm-sparse $head $covqual_fastgwa
     """
}
  process showFastGWAManhatten {
   memory params.other_process_memory
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


report_ch = report_fastlmm_ch.flatten().mix(pheno_report_ch.flatten())
                                     .mix(report_pca_ch.flatten())
				     .mix(report_plink_ch.flatten())
				     .mix(report_bolt_ch.flatten())
				     .mix(report_fastgwa_ch.flatten())
				     .mix(report_gemma_ch_GxE.flatten())
				     .mix(report_plink_gxe.flatten())
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





