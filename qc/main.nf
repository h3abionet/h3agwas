#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *      Jean-Tristan Brandenburg
 *      Shaun Aron
 *      Rob Clucas
 *      Eugene de Beste
 *      Scott Hazelhurst
 *      Anmol Kiran
 *      Lerato Magosi
 *      Abayomi Mosaku
 *
 *  On behalf of the H3ABionet Consortium
 *  2015-2020
 *
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 *
 *(C) University of the Witwatersrand, Johannesburg, 2016-2020
 *    on behalf of the H3ABioNet Consortium
 *This is licensed under the MIT Licence. See the "LICENSE" file for details
 */

//---- General definitions --------------------------------------------------//

import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;


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

    



def helps = [ 'help' : 'help' ]
params.batch = "0"
params.phenotype="0"
params.samplesheet   = "0"
if (params.idpat ==  "0") 
    idpat   = "(.*)"
else
    idpat   = params.idpat


params.help = false
K = "--keep-allele-order"

def params_help = new LinkedHashMap(helps)


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



nullfile = [false,"False","false", "FALSE",0,"","0","null",null]



def checkColumnHeader(fname, columns) {
  if (workflow.profile == "awsbatch") return;
  if (fname.toString().contains("s3://")) return;
  if (nullfile.contains(fname)) return;
  new File(fname).withReader { line = it.readLine().tokenize() }  
  problem = false; 
  columns.each { col -> 
    if (! line.contains(col) ) {
      println "The file <$fname> does not contain the column <$col>";
      problem=true;
    }
    if (problem)
      System.exit(2)
  }
}


def checkSampleSheet(fname)  {
  if (workflow.profile == "awsbatch") return;
  if (fname.contains("s3://") )return;
  if (nullfile.contains(fname) || fname.contains(".xls")) return;
  new File(fname).withReader { line = it.readLine()}  
  problem  = false
  prob_str = ""
  if (! line.contains(",")) {
    problem = true;
    prob_str = "If given as a CSV file, it must be comma-separated\n";
  }
  headers = line.tokenize(",")
  headers.each { println it}
  if (!(headers.contains("Institute Sample Label") || 
      (headers.contains("Sample Plate") && headers.contains("Well")))) {
    problem= true
    prob_str = prob_str + "Column headers must include 'Institute Sample Label'  or both 'Sample Plate' and 'Well'"
  }
  if (problem)  {
    println "There's a problem with the sample sheet <$fname>."
    println prob_str;
    System.exit(1)
  }
}        

if (nullfile.contains(params.samplesheet))
     samplesheet = "0"
   else {
      samplesheet = params.samplesheet
     checkSampleSheet(samplesheet)
   }

idfiles = [params.batch,params.phenotype]
idfiles.each { checkColumnHeader(it,['FID','IID']) }

println "The batch file is ${params.batch}"




nextflowversion =nextflow.version


if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"





// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}




params.queue    = 'batch'
params.remove_on_bp  = 1





max_plink_cores = params.max_plink_cores 
plink_mem_req   = params.plink_mem_req
other_mem_req   = params.other_mem_req
pi_hat          = params.pi_hat
super_pi_hat    = params.super_pi_hat
cut_diff_miss   = params.cut_diff_miss
f_lo_male       = params.f_lo_male
f_hi_female     = params.f_hi_female
remove_on_bp    = params.remove_on_bp

allowed_params= ["AMI","accessKey","batch","batch_col","bootStorageSize","case_control","case_control_col", "chipdescription", "cut_het_high","cut_get_low","cut_maf","cut_mind","cut_geno","cut_hwe","f_hi_female","f_lo_male","cut_diff_miss","cut_het_low", "help","input_dir","input_pat","instanceType","manifest", "maxInstances", "max_plink_cores","high_ld_regions_fname","other_mem_req","output", "output_align", "output_dir","phenotype","pheno_col","pi_hat", "plink_mem_req","region","reference","samplesheet", "scripts","secretKey","sexinfo_available", "sharedStorageMount","strandreport","work_dir","max_forks","big_time","super_pi_hat","samplesize","idpat","newpat","access-key","secret-key","instance-type","boot-storage-size","max-instances","shared-storage-mount","gemma_num_cores","remove_on_bp","queue","data","pheno","gc10"]



params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
    	println "Check $parm  ************** is it a valid parameter -- are you using one rather than two - signs or vice-versa";
  }
}

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


// This method first checks that the data file has the stated column 
// If so, it creates a channel for it
// NB: if the file is in S3 we cannot do the test since Groovy does not
// allow us to access the file directly
def getSubChannel = { parm, parm_name, col_name ->
  if (parm.toString().contains("s3://")) {
    println "The file <$parm> is in S3 so we cannot do a pre-check";
    return Channel.fromPath(parm);
  }
  if ((parm==0) || (parm=="0") || (parm==false) || (parm=="false")) {
    filename = "emptyZ0${parm_name}.txt";
    new File(filename).createNewFile()  
    new_ch = Channel.fromPath(filename);
  } else {
    if (! file(parm).exists()) {
     error("\n\nThe file <$parm> given for <params.${parm_name}> does not exist")
    } else {
      def line  
      new File(parm).withReader { line = it.readLine() }  
      fields = line.split()
      if (! fields.contains(col_name))
	  error("\n\nThe file <$parm> given for <params.${parm_name}> does not have a column <${col_name}>\n")
    }
    new_ch = Channel.fromPath(parm);
  }
  return new_ch;
}




if (params.case_control) {
  ccfile = params.case_control
  Channel.fromPath(ccfile).into { cc_ch; cc2_ch }
  col    = params.case_control_col
  diffpheno = "--pheno cc.phe --pheno-name $col"
  if (params.case_control.toString().contains("s3://")) {
       println "Case control file is in s3 so we can't check it"
  } else 
  if (! file(params.case_control).exists()) {
     error("\n\nThe file <${params.case_control}> given for <params.case_control> does not exist")
    } else {
      def line  
      new File(params.case_control).withReader { line = it.readLine() }  
      fields = line.split()
      if (! fields.contains(params.case_control_col))
	  error("\n\nThe file <${params.case_control}> given for <params.case_control> does not have a column <${params.case_control_col}>\n")
    }

} else {
  diffpheno = ""
  col = ""
  cc_ch  = Channel.value.into("none").into { cc_ch; cc2_ch }
}


phenotype_ch = getSubChannel(params.phenotype,"pheno",params.pheno_col)
batch_ch     = getSubChannel(params.batch,"batch",params.batch_col)




//---- Modification of variables for pipeline -------------------------------//


/* Define the command to add for plink depending on whether sexinfo is
 * available or not.
 */





if ( nullfile.contains(params.sexinfo_available) ) {
  sexinfo = "--allow-no-sex"
  extrasexinfo = ""
  println "Sexinfo not available, command --allow-no-sex\n"
} else {
  sexinfo = ""
  extrasexinfo = "--must-have-sex"
  println "Sexinfo available command"
}



/* Get the input files -- could be a glob
 * We match the bed, bim, fam file -- order determined lexicographically
 * not by order given, we check that they exist and then 
 * send the all the files to raw_ch and just the bim file to bim_ch */
inpat = "${params.input_dir}/${params.input_pat}"



if (inpat.contains("s3://")) {
  print "Here"
  this_checker = { it -> return it}
} else {
  this_checker = checker
}

Channel
   .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true)
    { file -> file.baseName }  \
      .ifEmpty { error "No matching plink files" }        \
      .map { a -> [this_checker(a[1]), this_checker(a[2]), this_checker(a[3])] }\
      .multiMap {  it ->
         raw_ch: it
         bim_ch: it[1] 
         inpmd5ch : it
       }.set {checked_input}

  


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



// Generate MD5 sums of output files
process inMD5 {
  input:
     file plink from checked_input.inpmd5ch
  output:
     file(out) into report_inpmd5_ch
  echo true
  script:
       bed = plink[0]
       bim = plink[1]
       fam = plink[2]
       out  = "${plink[0].baseName}.md5"
       template "md5.py"
}


println samplesheet
if (samplesheet != "0")  {
  sample_sheet_ch = file(samplesheet)

  process sampleSheet {
    input:
     file(sheet) from sample_sheet_ch
    output:
     file("poorgc10.lst") into poorgc10_ch
     file("plates") into report_poorgc10_ch
    script:
     """
       mkdir -p plates
       sampleqc.py $sheet ${params.gc10} "${idpat}"  poorgc10.lst plates/crgc10.tex
      """
    }
} else {

  process noSampleSheet {
    output:
     file("poorgc10.lst") into poorgc10_ch
     file("plates") into report_poorgc10_ch
    script:
     """
       mkdir -p plates
       sampleqc.py 0 0 0 poorgc10.lst plates/crgc10.tex
      """
    }

}  
  



//---- Start Pipeline -------------------------------------------------------//




/* Process to find duplicates. *
 * Inputs:
 * - bim: the bim file
 * Outputs:
 * - duplicates.snps    : A possibly empty file with a list of SNPs
 */
process getDuplicateMarkers {
  memory other_mem_req
  publishDir params.output_dir, pattern: "*dups", \
             overwrite:true, mode:'copy'
  input:
    file(inpfname) from checked_input.bim_ch
  output:
    file("${base}.dups") into duplicates_ch
  script:
     base     = inpfname.baseName
     outfname = "${base}.dups"
     template "dups.py"
}



/*  Process to remove duplicate SNPs.
 * Inputs
 *   -  raw files from from user-specified data
 *   -  list of duplicates comes from getDuplicateMarkers
 * Outputs:
 *   nodups.{bed,bim,fam} (PLINK file without duplicates) and
 *   qc.log log file
 */

process removeDuplicateSNPs {
  memory plink_mem_req
  input:
   tuple file(bed), file(bim), file(fam) from checked_input.raw_ch
   file(dups) from  duplicates_ch

  output:
    tuple  file("${nodup}.bed"),file("${nodup}.bim"),file("${nodup}.fam")\
    into (qc1_ch,qc1B_ch,qc1C_ch,qc1D_ch,qc1E_ch)
    tuple file("${base}.orig"), file(dups) into report_dups_ch
    file ("${nodup}.lmiss") into snp_miss_ch
    file ("${nodup}.imiss") into (ind_miss_ch1, ind_miss_ch2)
  script:
   base    = bed.baseName
   nodup   = "${base}-nd"
   """
    plink $K --bfile $base $sexinfo $extrasexinfo --exclude $dups --missing --make-bed --out $nodup
    wc -l ${base}.bim > ${base}.orig
    wc -l ${base}.fam >> ${base}.orig
   """
}


missingness = [0.01,0.03,0.05]  // this is used by one of the templates

if (extrasexinfo == "--must-have-sex") {


   /* Detailed analysis of X-chromosome */
   process getX {
     memory other_mem_req
     input:
       file(plink) from qc1D_ch
      output:
       file("X*") into X_chr_ch   
      script:
      base = plink[0].baseName
      """
	if [[ `grep  "^23" *bim`  ]];  then
	   plink --bfile $base --chr 23 --geno 0.04 --make-bed --out X
	else
	   echo ""
	   echo "----------------------------------------"
	   echo "There are no X-chromosome SNPs in this data "
	   echo "it does not make sense to check for sex"
	   echo "set sexinfo_available to false"
	   echo "----------------------------------------"
	   echo ""
	   exit 23 
	   touch X.bed X.bim X.fam EMPTYX
	fi
      """
   }



   process analyseX {
     memory other_mem_req
     input:
       file(xchr) from X_chr_ch
     output:
       file(out) into x_analy_res_ch // batchReport 
     script:
	x = xchr[0].baseName
	out = "x.pkl"
	template "xCheck.py"
   }
} else {


  x_analy_res_ch = Channel.fromPath("0")
   

}
/* Process to identify individual discordant sex information.
 * results are put in the output directory
 * Also does HWE
 */
process identifyIndivDiscSexinfo {
  memory plink_mem_req
  input:
     file(plinks) from qc1B_ch

  publishDir params.output_dir, overwrite:true, mode:'copy'

  output:
     file(logfile) into  (report_failed_sex_ch, failed_sex_ch1)
     tuple file(imiss), file(lmiss),file(sexcheck_report) into batchrep_missing_ch
     file("${base}.hwe") into hwe_stats_ch
  validExitStatus 0, 1
  script:
    base = plinks[0].baseName
    logfile= "${base}.badsex"
    sexcheck_report = "${base}.sexcheck"
    imiss  = "${base}.imiss"
    lmiss  = "${base}.lmiss"
    if (params.sexinfo_available == true)
    """
       plink $K --bfile $base --hardy --check-sex $f_hi_female $f_lo_male --missing  --out $base
       head -n 1 ${base}.sexcheck > $logfile
       grep  'PROBLEM' ${base}.sexcheck >> $logfile
    """
    else
     """
     plink --bfile $base  --hardy --missing  --out $base
     echo 'FID IID STATUS' > $sexcheck_report
     echo 'No sex information available to check'  > $logfile
     """
}



process generateSnpMissingnessPlot {
  memory other_mem_req
  input:
      file(lmissf) from snp_miss_ch
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
     file(output) into report_snpmiss_ch

  echo true
  script:
    input  = lmissf
    base   = lmissf.baseName
    label  = "SNPs"
    output = "${base}-snpmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}


process generateIndivMissingnessPlot {
  memory other_mem_req
  input:
      file(imissf) from ind_miss_ch1
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file(output) into report_indmisspdf_ch
  script:
    input  = imissf
    base   = imissf.baseName
    label  = "samples"
    output = "${base}-indmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}
 
process getInitMAF {
  memory plink_mem_req
  input:
     file(plink) from qc1C_ch
  output:
     file("${newbase}.frq") into init_freq_ch
  script:
    base = plink[0].baseName
    newbase = base.replace(".","_")
    """
    plink --bfile $base --freq --out $newbase
    """
}


process showInitMAF {
  memory other_mem_req
  input:
     file(freq) from init_freq_ch
  output:
     tuple file("${base}.pdf"), file("${base}.tex") into report_initmaf_ch
  script:
    base = freq.baseName+"-initmaf"
    base = base.replace(".","_")
    template "showmaf.py"
}

process showHWEStats {
  memory other_mem_req
  input:
     file(hwe) from hwe_stats_ch
  output:
     tuple file("${base}.pdf"), file("${base}-qq.pdf"), file("${base}.tex") into report_inithwe_ch
  script:
    base = hwe.baseName+"-inithwe"
    base = base.replace(".","_")
    template "showhwe.py"
}


process removeQCPhase1 {
  memory plink_mem_req
  input:
    tuple file(bed), file(bim), file(fam) from qc1_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file("${output}*.{bed,bim,fam}") into (qc2A_ch,qc2B_ch,qc2C_ch,qc2D_ch)
     tuple file("qc1.out"), file("${output}.irem") into report_qc1_ch
  script:
     base=bed.baseName
     output = "${base}-c".replace(".","_")
     """
     # remove really realy bad SNPs and really bad individuals
     plink $K --autosome --bfile $base $sexinfo --mind 0.1 --geno 0.1 --make-bed --out temp1
     plink $K --bfile temp1  $sexinfo --mind $params.cut_mind --make-bed --out temp2
     plink $K --bfile temp2  $sexinfo --geno $params.cut_geno --make-bed --out temp3
     plink $K --bfile temp3  $sexinfo --maf $params.cut_maf --make-bed --out temp4
     plink $K --bfile temp4  $sexinfo --hwe $params.cut_hwe --make-bed  --out $output 
     cat *log > logfile
     touch tmp.irem
     cat *.irem > ${output}.irem
     qc1logextract.py logfile ${output}.irem > qc1.out     
  """
}


// We do PCA on qc2 data because relatively few SNPs and individuals will be removed later and
// this is an expensive operation so we start early. Similarly for computing relatedness
process compPCA {
   cpus max_plink_cores
   memory plink_mem_req
   input:
      file plinks from qc2A_ch
   output:
      tuple file ("${prune}.eigenval"), file("${prune}.eigenvec") into (pcares, pcares1)
      tuple file ("${prune}.bed"), file("${prune}.bim"), file("${prune}.fam") into out_only_pcs_ch
   publishDir "${params.output_dir}/pca", overwrite:true, mode:'copy',pattern: "${prune}*"
   script:
      base = plinks[0].baseName
      prune= "${base}-prune".replace(".","_")
     """
     plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
     plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
     plink --bfile ${prune} --pca --out $prune 
     """
}

process drawPCA {
    memory other_mem_req
    input:
      tuple file(eigvals), file(eigvecs) from pcares
      file cc from cc2_ch
    output:
      tuple  file ("eigenvalue.pdf"), file(output) into report_pca_ch
    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.pdf"
    script:
      base=eigvals.baseName
      cc_fname = params.case_control
      // also relies on "col" defined above
      output="${base}-pca".replace(".","_")+".pdf"
      template "drawPCA.py"

}



if (params.high_ld_regions_fname != "")
   ldreg_ch=Channel.fromPath(params.high_ld_regions_fname)
   ////check
else
   ldreg_ch=Channel.value("dummy") //If not, create dummy channel






// Get which SNPs should be pruned for IBD
process pruneForIBD {
	// multi-threaded plink -- probably 2 core optimal, maybe 3
  cpus max_plink_cores
  memory plink_mem_req
  input:
    file plinks from qc2B_ch
    file ldreg  from ldreg_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file "${outf}.genome" into (find_rel_ch,batch_rel_ch)
  script:
    base   = plinks[0].baseName
    outf   =  base.replace(".","_")
    if (params.high_ld_regions_fname != "")
      range = " --exclude range $ldreg"
    else
      range =""
    """
     plink --bfile $base --threads $max_plink_cores --autosome $sexinfo $range --indep-pairwise 60 5 0.2 --out ibd
     plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --out ibd_prune
     plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --min $pi_hat --out $outf
      echo DONE
     """
}



// run script to find a tuple of individuals we can remove to ensure no relatedness
//  Future - perhaps replaced with Primus
process findRelatedIndiv {
  memory other_mem_req
  input:
     file (missing) from ind_miss_ch2
     file (ibd_genome) from find_rel_ch
  output:
     file(outfname) into (related_indivs_ch1,related_indivs_ch2, report_related_ch) 
  publishDir params.output_dir, overwrite:true, mode:'copy'
  script:
     base = missing.baseName
     outfname = "${base}-fail_IBD".replace(".","_")+".txt"
     template "removeRelInds.py"
}



process calculateSampleHeterozygosity {
   memory plink_mem_req
   input:
      file(nodups) from qc2C_ch

   publishDir params.output_dir, overwrite:true, mode:'copy'
   output:
      tuple file("${hetf}.het"), file("${hetf}.imiss") into (hetero_check_ch, plot1_het_ch)
      file("${hetf}.imiss") into missing_stats_ch
   script:
      base = nodups[0].baseName
      hetf = "${base}".replace(".","_")
   """
     plink --bfile $base  $sexinfo --het --missing  --out $hetf
   """
}



process generateMissHetPlot {
  memory other_mem_req
  input:
    tuple file(het), file(imiss) from plot1_het_ch
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file(output) into report_misshet_ch
  script:
    base = imiss.baseName
    output  = "${base}-imiss-vs-het".replace(".","_")+".pdf"
    template "missHetPlot.py"
}



// Find those who have bad heterozygosity
process getBadIndivsMissingHet {
  memory other_mem_req
  input:
    tuple file(het), file(imiss) from hetero_check_ch
  output:
    file(outfname) into (failed_miss_het, report_misshetremf_ch)
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.txt"
  script:
    base = het.baseName
    outfname = "${base}-fail_het".replace(".","_")+".txt"
    template "select_miss_het_qcplink.py"
}




process removeQCIndivs {
  memory plink_mem_req
  input:
    file(f_miss_het)     from failed_miss_het
    file(rel_indivs)     from related_indivs_ch1
    file (f_sex_check_f) from failed_sex_ch1
    file (poorgc)        from poorgc10_ch
    tuple file(bed), file(bim), file(fam) from qc2D_ch
  output:
     file("${out}.{bed,bim,fam}") into\
        (qc3A_ch, qc3B_ch)
  script:
   base = bed.baseName
   out  = "${base}-c".replace(".","_")
    """
     cat $f_sex_check_f $rel_indivs $f_miss_het $poorgc | sort -k1 | uniq > failed_inds
     plink $K --bfile $base $sexinfo --remove failed_inds --make-bed --out $out
     mv failed_inds ${out}.irem
  """
}


mperm_header=" CHR                               SNP         EMP1         EMP2 "

// Find differential missingness between cases and controls; also compute HWE scores
process calculateSnpSkewStatus {
  memory plink_mem_req
  cpus max_plink_cores
  input:
    file(plinks) from qc3A_ch.combine(cc_ch)
  output:
    file "${base}.missing" into clean_diff_miss_plot_ch1
    file mperm into clean_diff_miss_ch2
    file "${base}.hwe" into hwe_scores_ch
  script:
   base  = plinks[0].baseName
   out   = base.replace(".","_")
   mperm = "${base}.missing.mperm"
   phe   = plinks[3]
   """
    cp $phe cc.phe
    plink --threads ${max_plink_cores} --autosome --bfile $base $sexinfo $diffpheno --test-missing mperm=10000 --hardy --out $out
    if ! [ -e $mperm ]; then
       echo "$mperm_header" > $mperm
    fi

   """
}


process generateDifferentialMissingnessPlot {
   memory other_mem_req
   input:
     file clean_missing from clean_diff_miss_plot_ch1
   publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
   output:
      file output into report_diffmissP_ch
   script:
       input = clean_missing
       base  = clean_missing.baseName.replace(".","_").replace("-nd","")
       output= "${base}-diff-snpmiss_plot.pdf"
       template "diffMiss.py"

 }


// Find those SNPs that have diff missingness in cases & controls
process findSnpExtremeDifferentialMissingness {
  memory other_mem_req
  input:
    file clean_missing from clean_diff_miss_ch2
  echo true
  output:
     tuple val(base), file(failed) into bad_snps_ch
     file(failed) into report_diffmiss_ch
     file(failed) into skewsnps_ch
  script:
    cut_diff_miss=params.cut_diff_miss
    missing = clean_missing
    base     = missing.baseName.replace("-.*","").replace(".","_")
    probcol = 'EMP2'  // need to change if we don't use mperm
    failed   = "${base}-failed_diffmiss.snps"
    template "select_diffmiss_qcplink.py"
}


process removeSkewSnps {
  memory plink_mem_req
  input:
    file (plinks) from qc3B_ch
    file(failed) from skewsnps_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'  
  output:
    tuple file("${output}.bed"), file("${output}.bim"), file("${output}.fam"), file("${output}.log") \
      into (qc4A_ch, qc4B_ch, qc4C_ch,  report_cleaned_ch)
  script:
  base = plinks[0].baseName
  output = params.output.replace(".","_")
  """
  plink $K --bfile $base $sexinfo --exclude $failed --make-bed --out $output
  """
}




/*
process convertInVcf {
   memory plink_mem_req
   cpus max_plink_cores
   input :
     tuple file(bed), file(bim), file(fam), file (log) from qc4A_ch
   publishDir params.output_dir, overwrite:true, mode:'copy'
   output :
    file("${base}.vcf")  
   script:
     base= bed.baseName
     """
     plink --bfile ${base} --threads ${max_plink_cores} --recode vcf --out $base
     """

}*/


process calculateMaf {
  memory plink_mem_req
  input:
    tuple  file(bed), file(bim), file(fam), file(log) from qc4C_ch

  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.frq"

  output:
    file "${base}.frq" into maf_plot_ch

  script:
    base = bed.baseName
    out  = base.replace(".","_")
    """
      plink --bfile $base $sexinfo  --freq --out $out
    """
}



process generateMafPlot {
  memory other_mem_req
  input:
    file input from maf_plot_ch
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file(output) into report_mafpdf_ch

  script:
    base    = input.baseName
    output  = "${base}-maf_plot.pdf"
    template "mafplot.py"
}




// Find HWE scores of each SNP
process findHWEofSNPs {
  memory other_mem_req
  input:
     file hwe from hwe_scores_ch
  output:
     file output  into unaff_hwe

  script:
    base   = hwe.baseName.replace(".","_")
    output = "${base}-unaff.hwe"
    """
      head -1 $hwe > $output
      grep 'UNAFF' $hwe >> $output
    """
}

process generateHwePlot {
  memory other_mem_req
  input:
    file unaff from unaff_hwe
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file output into report_hwepdf_ch

  script:
    input  = unaff
    base   = unaff.baseName.replace(".","_")
    output = "${base}-hwe_plot.pdf"
    template "hweplot.py"
}




// Generate MD5 sums of output files
process outMD5 {

  input:
     tuple file(bed), file(bim), file(fam), file(log) from qc4B_ch
  output:
     file(out) into report_outmd5_ch
  echo true
  script:
       out  = "${bed.baseName}.md5"
       template "md5.py"
}





process batchProc {
  memory plink_mem_req
  input:
    tuple file(eigenval), file(eigenvec) from pcares1
    tuple file(imiss), file(lmiss), file(sexcheck_report) from batchrep_missing_ch
    file "pheno.phe" from phenotype_ch    // staged input file
    file "batch.phe" from batch_ch        // staged input file
    file genome    from batch_rel_ch    // pruneForIBD
    file pkl       from x_analy_res_ch  // analyseX
    file rem_indivs from related_indivs_ch2 // findRel
  publishDir params.output_dir, pattern: "*{csv,pdf}", \
             overwrite:true, mode:'copy'
  output:
      file("${base}-batch.tex")      into report_batch_report_ch
      tuple file("*.csv"), file("*pdf") into report_batch_aux_ch // need to stage
  script:
    phenotype = "pheno.phe"
    batch = "batch.phe"
    base = eigenval.baseName
    batch_col = params.batch_col
    pheno_col = params.pheno_col
    template "batchReport.py"
}




repnames = ["dups","cleaned","misshet","mafpdf","snpmiss","indmisspdf","failedsex","misshetremf","diffmissP","diffmiss","pca","hwepdf","related","inpmd5","outmd5","batch"]



process produceReports {
  memory other_mem_req
  label 'latex'
  input:
    tuple file(orig), file (dupf) from report_dups_ch
    tuple file(cbed), file(cbim), file(cfam), file(ilog) from report_cleaned_ch
    file(missingvhetpdf) from report_misshet_ch
    file(mafpdf)         from report_mafpdf_ch
    file(snpmisspdf)     from report_snpmiss_ch
    file(indmisspdf)     from report_indmisspdf_ch
    file(fsex)           from report_failed_sex_ch
    file(misshetremf)    from report_misshetremf_ch
    file(diffmisspdf)    from report_diffmissP_ch
    file(diffmiss)       from report_diffmiss_ch
    tuple file(eigenvalpdf),file(pcapdf)         from report_pca_ch
    file(hwepdf)         from report_hwepdf_ch
    file(rel_indivs)     from report_related_ch
    file(inpmd5)         from report_inpmd5_ch
    file(outmd5)         from report_outmd5_ch
    tuple file(initmafpdf), file(initmaftex) from report_initmaf_ch
    tuple file(inithwepdf), file(inithweqqpdf), file(inithwetex) from report_inithwe_ch
    tuple file(qc1), file(irem)  from report_qc1_ch
    file(batch_tex)  from report_batch_report_ch
    file(poorgc)     from report_poorgc10_ch
    tuple file(bpdfs), file(bcsvs) from report_batch_aux_ch
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file("${base}.pdf") into final_ch
   script:
     base = params.output
     config_text = getConfig()
     template "qcreport.py"
}



final_ch.subscribe { b=it.baseName; 
  new File("emptyZ0batch.txt").delete();
  new File("emptyZ0pheno.txt").delete();
  new File("xxemptyZ0pheno.txt").delete();
  println "The output report is called ${params.output_dir}/${b}.pdf"}
