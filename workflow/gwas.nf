#!/usr/bin/env nextflow

/*
 * Authors       : 
 *
 *      
 *      Shaun Aron
 *   	Rob Clucas
 *      Eugene de Beste
 *      Scott Hazelhurst
 *      Abayomi Mosaku
 *      Anmol Kiran
 *      Lerato Magosi
 *      
 *  On behalf of the H3ABionet Consortium
 *  2015-2016
 * 
 *
 * Description  : Nextflow pipeline for Wits GWAS.
 * 
 */

//---- General definitions --------------------------------------------------//



/* Defines the directory where the plink 2 input binary files are. 
 *
 * NOTE: This must be a relative path, from where the pipeline is run.
 * and should end with a slash
 */
import java.nio.file.Paths

def helps = [ 'help' : 'help' ]


def params_help = new LinkedHashMap(helps)


params.work_dir   = "$HOME/h3agwas"
params.input_dir  = "${params.work_dir}/input"  
params.output_dir = "${params.work_dir}/output"



/* Defines the path where any scripts to be executed can be found.
 */

params.scripts   = "${params.work_dir}/scripts"

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

//---- Cutoff definitions ---------------------------------------------------//

/* Defines the cutoffs for the heterozygosity. Standard cutoff +- 3sd from 
 * mean)
 */
params.cut_het_high = 0.343
params.cut_het_low  = 0.254

/* Defines the cutoff for missingness. Using standard cutoff -- 3 - 7%.
 */
params.cut_miss      = 0.05
params.cut_diff_miss = 0.05;


/* Defines the cutoff for the SNP minor allele frequency.
 */
params.cut_maf        = 0.01

/* Defines the cutoff for SNP missingness.
 */
params.cut_mind     = 0.01
params.cut_geno     = 0.01

/* Defines the cutoff for the SNP Hardy Weinburg deviation.
 */
params.cut_hwe        = 0.008

params.help = false
    /* cut-off for relatedness */

params.pi_hat = 0.04
pi_hat=params.pi_hat

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

// TODO  : fix as nice Groovy
bim = Paths.get(params.input_path,"${params.plink_fname}.bim").toString()



// Prepends scripts directory path to the argument given
def path = { 
   fn ->
     Paths.get(params.scripts,fn).toString()
}


// Checks if the file exists
def checker = { fn -> 
   if (fn.exists()) 
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
  }


//------------

// Creating two channels with the file names and at the same time
// checking file existence



bed = Paths.get(params.input_path,"${params.data_name}.bed").toString()
println "Asking for $bed"
bim = Paths.get(params.input_path,"${params.data_name}.bim").toString()
fam = Paths.get(params.input_path,"${params.data_name}.fam").toString()
bim_ch = Channel.fromPath(bim).map checker
Channel
    .from(file(bed),file(bim),file(fam))
    .buffer(size:3)
    .map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }
    .set { raw_ch }


//---- Start Pipeline -------------------------------------------------------//

/* Process to find duplicates. * 
 * Inputs:
 * - bim: the bim file
 * Outputs:
 * - duplicates.snps    : A possibly empty file with a list of SNPs
 */
process getDuplicateMarkers { 
  input:
    set file("raw.bim") from bim_ch
  output:
    set file("duplicates.snps") into remove_ch
  script:
    dup_find = path("dups.py")
    """
      $dup_find raw.bim duplicates.snps
    """
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
  input:
    set file(bed), file(bim), file(fam) from raw_ch
    set file('duplicates.snps') from remove_ch

  publishDir params.publish, pattern: "0002-dups.log", overwrite:true, mode:'copy'

  output:
    set  file('nodups.bed'),file('nodups.bim'),file('nodups.fam') into ready_ch
    file ('0002-dups.log')
  script:
   base=bed.baseName
   """
    plink --bfile ${base} $sexinfo --exclude duplicates.snps --make-bed --out nodups >> qc.log
    mv nodups.log 0002-dups.log
   """
}


// Now make copies of  the ready_ch so work can happen in parallel

(sex_check_ch,missing_ch,het_ch,ibd_prune_ch,remove_inds_ch) = [1,2,3,4,5].collect {Channel.create()}
ready_ch.into(sex_check_ch,missing_ch,het_ch,remove_inds_ch,ibd_prune_ch) 



/* Process to identify individual discordant sex information.
 * results are put in the output directory
 */
process identifyIndivDiscSexinfo {
  input:
     set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from sex_check_ch
  
  publishDir params.publish, overwrite:true, mode:'copy'

  output:
     file '0010-failed.sexcheck' into failed_sex_check
  script:
  if (params.sexinfo_available == "true")
  """
       plink --bfile nodups --check-sex  --out nodups
       if grep -Rn 'PROBLEM' nodups.sexcheck > 0010-failed.sexcheck; then
         echo 'Discordant sex info found'
       else                                                      
         echo 'No discordant sex info found'
       fi

  """
  else
    "echo 'No sex information available to check'  > 0010-failed.sexcheck"

}


// Find missingness statistics for the plink file
process calculateSampleMissing {
  input:
     set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from missing_ch

  publishDir params.publish, overwrite:true, mode:'copy'

  output:
     file("0020.imiss") into calc_missing_ch
  """
    plink --bfile nodups $sexinfo --missing --out 0020
  """
}


// Will use the missingness statictis twice -- once for a report
// once for removing data
plot1_ch_miss = Channel.create()
miss_het_ch   = Channel.create()
missing2_ch   = Channel.create()

calc_missing_ch.into(plot1_ch_miss,missing2_ch,miss_het_ch)



    // But for the  moment let's deal with heterozygosity

process calculateSampleHetrozygosity {
   input:
      set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from het_ch

   publishDir params.publish, overwrite:true, mode:'copy'

   output:
      file("0030.het") into sample_hetero_result
   script:
   """ 
     plink --bfile nodups $sexinfo --het  --out 0030
   """
}

channels = (hetero_check_ch, plot1_ch_het)= [1,2,3,4].collect {Channel.create()}
sample_hetero_result.into(hetero_check_ch,plot1_ch_het)

process generateMissHetPlot {
  errorStrategy 'ignore'

  input:
  file 'qcplink.imiss' from plot1_ch_miss
  file 'qcplink.het'   from plot1_ch_het   

  publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"

  output:
    file('*.pdf')   into pictures_ch

  script:
    plotscript = path("miss_het_plot_qcplink.R")
    """
     $plotscript qcplink.imiss qcplink.het pairs.imiss-vs-het.pdf meanhet_plot.pdf
    """
}


// Find those who have too high missingness, or bad heterozygosity
process getBadIndivs_Missing_Het {
  errorStrategy 'ignore'

  input:
  file 'qcplink.imiss' from miss_het_ch
  file 'qcplink.het'   from hetero_check_ch

  output:
    file('fail_miss_het_qcplink.txt') into failed_miss_het

  script:
    selectscript = path("select_miss_het_qcplink.pl")
    """
     $selectscript $params.cut_het_high $params.cut_het_low $params.cut_miss \
                     qcplink.imiss qcplink.het fail_miss_het_qcplink.txt
    """
}

/* We are going to check for related individuals and remove them */

// first, in computing relatedness do we ignore high LD regions? 

if (params.high_ld_regions_fname != "")
   ldreg_ch=Channel.fromPath(params.plink_inputpath+params.high_ld_regions_fname)
else 
   ldreg_ch=Channel.value("dummy") //If not, create dummy channel


// Get which SNPs should be pruned for IBD
process pruneForIBD {
	// multi-threaded plink -- probably 2 core optimal, maybe 3
  input:
    set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from ibd_prune_ch
    file ldreg    from ldreg_ch
  output:
  //set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') into ibd
    file 'ibd_min_thresh.genome' into sort_ibd_ch1,sort_ibd_ch2
  script:
    if (params.high_ld_regions_fname != "")
      range = "--range --exclude $ldreg"
    else
      range =""
    """
      plink --bfile nodups --autosome $sexinfo $range --indep-pairwise 50 5 0.2 --out ibd
      plink --bfile nodups --autosome $sexinfo --extract ibd.prune.in --genome --out ibd_prune
      plink --bfile nodups --autosome $sexinfo --extract ibd.prune.in --genome --min $pi_hat --out ibd_min_thresh
      echo DONE
     """

}

//(sort_ibd_ch1,sort_ibd_ch2) = [Channel.create(),Channel.create()]
//ibd_min_genome.into(sort_ibd_ch1,sort_ibd_ch2)

    // HUH?
process sortByPiHat {
  input:
     file(ibd_min_genome) from sort_ibd_ch1
  output:
     file 'qcplink_ibd_min_thresh_sorted_pihat.txt'

  """
  sort -k10n ${ibd_min_genome} > qcplink_ibd_min_thresh_sorted_pihat.txt
  """
}

// run script to find related individuals
//  Future - perhaps replaced with Primus
process findRelatedIndiv {
  errorStrategy 'ignore'

  input:
     file missing from missing2_ch
     file ibd_genome from sort_ibd_ch2

  output:
     file 'fail_IBD_qcplink.txt' into related_indivs

  script:
    ibdscript = path("run_IBD_QC_qcplink.pl")
  """
    $ibdscript $missing $ibd_genome fail_IBD_qcplink.txt
  """
}


process removeQCIndivs {
  input:
    file failed_miss_het
    file failed_sexcheck_f from failed_sex_check 
    file related_indivs
    set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from remove_inds_ch
  output:
     set file("clean00.bed"),file("clean00.bim"),file("clean00.fam") into \
         (clean00_ch1,clean00_ch2,clean00_ch3, clean00_ch4)

  script:
  """
  cat $failed_sexcheck_f $related_indivs $failed_miss_het | sort -k1 | uniq > qcplink_failed_inds
  plink --bfile nodups $sexinfo --remove qcplink_failed_inds --make-bed --out clean00
  """
}



process calculateMaf {
  input:
    set file("clean00.bed"),file("clean00.bim"),file("clean00.fam") from clean00_ch1

  publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.frq"

  output:
    file 'clean00.frq' into maf_plot_ch

  script:
  """
    plink --bfile clean00 $sexinfo  --freq --out clean00
  """
}


process generateMafPlot {
  input:
  file 'clean00.frq' from maf_plot_ch

  publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"

  output:
    file 'maf_plot.pdf'

  script:
    plotscript = path("maf_plot_qcplink.R")
    """
      $plotscript clean00.frq maf_plot.pdf
    """
}



// Repeat computation of missingness on QCd data
process calculateSnpMissigness {
  input:
   set file('clean00.bed'),file('clean00.bim'),file('clean00.fam')  from clean00_ch2

  output:
   file 'clean00.lmiss' into clean_miss_plot_ch

  script:
  """
   plink --bfile clean00 $sexinfo --missing --out clean00
  """
}


process generateSnpMissingnessPlot {
  input:
    file 'clean00.lmiss' from clean_miss_plot_ch

  publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"

  output:
    file 'snpmiss_plot.pdf'

  script:
    plotscript = path("snpmiss_plot_qcplink.R")
    """
      $plotscript clean00.lmiss snpmiss_plot.pdf
    """
}

// Find differential missingness between cases and controls; also compute HWE scores
process calculateSnpSkewStatus {
  input:
    set file('clean00.bed'),file('clean00.bim'),file('clean00.fam')  from clean00_ch3
  output:
    file 'clean00.missing' into src_clean_diff_miss_ch
    file 'clean00.hwe' into hwe_scores_ch
  script:
   """
    plink --bfile clean00 $sexinfo --test-missing --hardy --out clean00
   """
}  

(clean_diff_miss_plot_ch1,clean_diff_miss_ch2)=[1,2].collect{Channel.create()}
src_clean_diff_miss_ch.into(clean_diff_miss_plot_ch1,clean_diff_miss_ch2)

process generateDifferentialMissingnessPlot {
   input:
     file "clean00.missing" from clean_diff_miss_plot_ch1

   publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"
   output:
      file 'snpmiss_plot.pdf' into snpmiss_plot_ch

   script:
   plotscript = path("diffmiss_plot_qcplink.R")
   """
      $plotscript clean00.missing snpmiss_plot.pdf
   """
 }


// Find those SNPs that have diff missingness in cases & controls
process findSnpExtremeDifferentialMissingness {
  input:
    file "clean00.missing" from clean_diff_miss_ch2
  output:
     file 'failed_diffmiss.snps' into bad_snps_ch
  script:
    cut_diff_miss=params.cut_diff_miss
    diffscript = path("select_diffmiss_qcplink.pl")
    """ 
     $diffscript $cut_diff_miss clean00.missing failed_diffmiss.snps
    """
}

// Find HWE scores of each SNP
process findHWEofSNPs {
  input:
     file 'clean00.hwe' from hwe_scores_ch
  output:
     file 'unaff.hwe'   into unaff_hwe

  script:
  """
   head -1 clean00.hwe > unaff.hwe
   grep 'UNAFF' clean00.hwe >> unaff.hwe
  """
}

process generateHwePlot {
  input:
    file 'unaff.hwe' from unaff_hwe
  publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file 'hwe_plot.pdf'

  script:
    plotscript = path("hwe_plot_qcplink.R")
    """
     $plotscript unaff.hwe hwe_plot.pdf
    """

}


process removeQCPhase1 {

  input:
    set file('clean00.bed'),file('clean00.bim'),file('clean00.fam')  from clean00_ch4
    file 'failed_diffmiss.snps' from bad_snps_ch
  publishDir params.publish, overwrite:true, mode:'copy'
  output:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam')  into clean01_ch

  script:
  """
  plink --bfile clean00 $sexinfo --exclude failed_diffmiss.snps --mind 0.2 --make-bed --out temp1
  plink --bfile temp1 $sexinfo --geno 0.2 --make-bed --out temp2
  plink --bfile temp2  $sexinfo \
        --autosome \
        --maf $params.cut_maf --mind $params.cut_mind --geno $params.cut_geno --hwe $params.cut_hwe \
         --make-bed --out cleaned
  """
}

process computePhase0 {

  input:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam')  from clean01_ch

  publishDir params.publish, overwrite:true, mode:'copy'
  output:
    set file('cleaned.*')  into comp_phase1_ch

  script:
  """
     plink --bfile cleaned --pca --assoc --adjust --out cleaned
  """
}


pictures_ch.subscribe { println "Drawn $it" }
comp_phase1_ch.subscribe { println "Done!!!" }




