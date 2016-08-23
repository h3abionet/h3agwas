#!/usr/bin/env nextflow

/*
 * Author       : Rob Clucas
 * Description  : Nextflow pipeline for Wits GWAS.
 * Modified by Scott Hazelhurst
 */

//---- General definitions --------------------------------------------------//


/* Defines the name of the mountpoint of the data directories in the docker
 * container. This is so that any scripts which run in the container and 
 * might need this info can run succesfully, and the user can specify the 
 * directory to each of the scripts.
 *
 * NOTE: The mountpoint is mounted in the container from the root directory,
 *       so specifying 'util' as the mount point mounts the data at '/util' in
 *       the container.
 */
params.dock_mpoint      = 'util'

/* Defines the directory where the plink 2 input binary files are. 
 *
 * NOTE: This must be a relative path, from where the pipeline is run.
 * and should end with a slash
 */
params.plink_inputpath  = "/$HOME/witsGWAS/gwasdata/"  
params.publish          = "/$HOME/witsGWAS/output"

/* Defines the path where any scripts to be executed can be found.
 *
 * NOTE: This must be a ralative path, from where the pipeline is run.
 */
/* example /*
params.script_path      = 'scripts'

/* Defines the names of the plink binary files in the plink directory 
 * (.fam, .bed, .bed).
 *
 * NOTE: This must be without the extension (so if A.fam, A.bed, ... 
 *       then use 'A').
 */
params.plink_fname      = 'raw-GWA-data'

/* Defines the name of the file with high LD region information.
 * 
 * NOTE: This can have/cannot have the extension, but should be in the 
 *       plink_inputpath specified above.
 */
params.high_ld_regions_fname = 'high_LD_regions.txt'

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




//---- Modification of variables for pipeline -------------------------------//

/* Define the command to add for plink depending on whether sexinfo is
 * available or not. Command is:
 * 
 * - No sexinfo availabele  : "--allow-no-sexinfo"   
 * - Sexinfo available      : ""
 */
if ( params.sexinfo_available == "false" ) {
  sexinfo = "--allow-no-sex"
  println "Sexinfo not available, command --allow-no-sex\n"
} else {
  sexinfo = ""
  println "Sexinfo available command"
}

import java.nio.file.Paths

println params.plink_inputpath

bed = Paths.get(params.plink_inputpath,"${params.plink_fname}.bed").toString()
bim = Paths.get(params.plink_inputpath,"${params.plink_fname}.bim").toString()
fam = Paths.get(params.plink_inputpath,"${params.plink_fname}.fam").toString()


def script = { 
   fn ->
         Paths.get(params.script_path,fn).toString()
}


def checker = { fn -> 
   if (fn.exists()) 
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
  }



bim_ch = Channel.fromPath(bim).map checker
raw_ch = Channel.from(file(bed),file(bim),file(fam)).buffer(size:3).map { a -> [checker(a[0]), checker(a[1]), checker(a[2])] }



//---- Start Pipeline -------------------------------------------------------//

/* Process to check for duplicates. * 
 * Inputs:
 * - bim: the bim file
 * Outputs:
 * - duplicats.snps    : A possibly empty file with a list of SNPs
 */
process getDuplicateMarkers { 
  input:
  set file("raw.bim") from bim_ch

  output:
  set file("duplicates.snps") into remove_ch

  script:
  dup_find = script("dups.py")
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
process removeDuplicateMarkers {
  input:
    set file(bed), file(bim), file(fam) from raw_ch
    set file('duplicates.snps') from remove_ch

  publishDir params.publish, pattern: "nodups.log", overwrite:true, mode:'move'

  output:
    set  file('nodups.bed'),file('nodups.bim'),file('nodups.fam') into ready_ch
    file ('0002.log')
  script:
   base=bed.baseName
   """
    plink --bfile ${base} $sexinfo --exclude duplicates.snps --make-bed --out nodups >> qc.log
    mv nodups.log 0002.log
   """
}


// Now make copies of  the ready_ch so work can happen in parallel

sex_check_ch  = Channel.create()
missing_ch    = Channel.create()
het_ch        = Channel.create()
ibd_prune_ch  = Channel.create()
remove_indivs_ch = Channel.create()
ready_ch.into(sex_check_ch, missing_ch, het_ch, ibd_prune_ch, remove_indivs_ch) 



/* Process to identify individual discordant sex information.
 * results are put in the output directory
 */
/*Check this process*/
process identifyIndivDiscSexinfo {
  input:
     set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from sex_check_ch
  
  publishDir params.publish, overwrite:true, mode:'move'

  output:
     file("0010.sexcheck") 
     file 'failed.sex' into failed_sex_check
  script:
  """
  if [[ ${params.sexinfo_available} == 'true' ]]; then 
       plink --bfile nodups --check-sex  --out 0010 
       if grep -Rn 'PROBLEM' 0010.sexcheck > failed.sex; then
         echo 'Discordant sex info found'
       else                                                      
         echo 'No discordant sex info found'
        fi

  else
       touch 0010.sexcheck
       echo "No sex information available to check"  > failed.sex
  fi
  """
}


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

plot1_ch_miss = Channel.create()
missing2_ch   = Channel.create()

calc_missing_ch.into(plot1_ch_miss,missing2_ch) 

process calculateSampleHetrozygosity {
   input:
      set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from het_ch

   publishDir params.publish, overwrite:true, mode:'copy'

   output:
      file("0030.het") into plot1_ch_het

   script:
   """ 
     plink --bfile nodups $sexinfo --het  --out 0030
   """
}


process generateMissHetPlot {
  errorStrategy 'ignore'

  input:
  file 'qcplink.imiss' from plot1_ch_miss
  file 'qcplink.het'   from plot1_ch_het   

  publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"

  output:
    file('*.pdf')   into pictures_ch
    file('fail_miss_het_qcplink.txt') into failed_miss_het

  script:
    plotscript = script("miss_het_plot_qcplink.R")
    selectscript = script("select_miss_het_qcplink.pl")
  println plotscript
  """
  $plotscript qcplink.imiss qcplink.het pairs.imiss-vs-het.pdf meanhet_plot.pdf
  $selectscript $params.cut_het_high $params.cut_het_low $params.cut_miss \
                     qcplink.imiss qcplink.het fail_miss_het_qcplink.txt
  """
}




ldreg_ch=Channel.fromPath(params.plink_inputpath+params.high_ld_regions_fname)

// Get which SNPs should be pruned for IBD
process pruneForIBD {
  input:
    set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from ibd_prune_ch
    file ldreg    from ldreg_ch
  output:
    set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') into ibd
    file('ibd_min_0_04.genome') into ibd_min_0_04_genome
  script:
    """
      plink --bfile nodups --autosome $sexinfo --exclude $ldreg --range --indep-pairwise 50 5 0.2 --out ibd
      plink --bfile nodups --autosome $sexinfo --extract ibd.prune.in --genome --out ibd_prune
      plink --bfile nodups --autosome $sexinfo --extract ibd.prune.in --genome --min 0.04 --out ibd_min_0_04
     """
}

(sort_ibd_ch1,sort_ibd_ch2) = [Channel.create(),Channel.create()]
ibd_min_0_04_genome.into(sort_ibd_ch1,sort_ibd_ch2)

process sortByPiHat {
  input:
     file('ibd_min_0_04.genome') from sort_ibd_ch1
  output:
     file 'qcplink_ibd_min_0_04_sorted_pihat.txt'

  """
  sort -k10n ibd_min_0_04.genome > qcplink_ibd_min_0_04_sorted_pihat.txt
  """
}

process filterRelatedIndiv {
  errorStrategy 'ignore'

  input:
     file missing from missing2_ch
     file ibd_genome from sort_ibd_ch2

  output:
     file 'fail_IBD_qcplink.txt'

  script:
    ibdscript = script("run_IBD_QC_qcplink.pl")
  """
    $ibdscript $missing $ibd_genome fail_IBD_qcplink.txt
  """
}


process removeQCIndivs {
  input:
    file failed_miss_het
    file failed_sexcheck_f from failed_sex_check 
    set file('nodups.bed'),file('nodups.bim'),file('nodups.fam') from remove_indivs_ch
  output:
     set file("clean00.bed"),file("clean00.bim"),file("clean00.fam") into clean00_ch

  script:
  """
  ls > allfiles
  cat $failed_sexcheck_f $failed_miss_het | sort -k1 | uniq > qcplink_failed_inds
  plink --bfile nodups $sexinfo --remove qcplink_failed_inds --make-bed --out clean00
  """
}

(clean00_ch1, clean00_ch2, clean00_ch3, clean00_ch4)  = [1,2,3,4].collect {Channel.create()}
clean00_ch.into(clean00_ch1,clean00_ch2,clean00_ch3, clean00_ch4)

process calculateMaf {
  input:
     set file("clean00.bed"),file("clean00.bim"),file("clean00.fam") from clean00_ch1

  publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"

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
    plotscript = script("maf_plot_qcplink.R")
  """
    $plotscript clean00.frq maf_plot.pdf
  """
}



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
    plotscript = script("snpmiss_plot_qcplink.R")
  """
    $plotscript clean00.lmiss snpmiss_plot.pdf
  """
}

// Find differential missingness and HWE scores
process calculateSnpSkewStatus {
  input:
   set file('clean00.bed'),file('clean00.bim'),file('clean00.fam')  from clean00_ch3
  output:
   file 'clean00.missing' into src_clean_diff_miss_plot_ch
   file 'clean00.hwe' into hwe_scores_ch
  script:
  """
    plink --bfile clean00 $sexinfo --test-missing --hardy --out clean00
  """
}  

(clean_diff_miss_plot_ch1,clean_diff_miss_plot_ch2)=[1,2].collect{Channel.create()}
src_clean_diff_miss_plot_ch.into(clean_diff_miss_plot_ch1,clean_diff_miss_plot_ch2)

process generateDifferentialMissingnessPlot {
   input:
     file "clean00.missing" from clean_diff_miss_plot_ch1

    publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"
    output:
      file 'snpmiss_plot.pdf' into snpmiss_plot_ch

    script:
    plotscript = script("diffmiss_plot_qcplink.R")
     """
      $plotscript clean00.missing snpmiss_plot.pdf
     """
 }



process findSnpExtremeDifferentialMissingness {
  input:
    file "clean00.missing" from clean_diff_miss_plot_ch2
  output:
     file 'failed_diffmiss.snps' into bad_snps_ch
  script:
    cut_diff_miss=params.cut_diff_miss
    diffscript = script("select_diffmiss_qcplink.pl")
    """ 
     $diffscript $cut_diff_miss clean00.missing failed_diffmiss.snps
    """
}

process findSnpsExtremeHweDeviations {
  input:
     file 'clean00.hwe' from hwe_scores_ch
  output:
     file 'unaff.hwe'   into unaff

  script:
  """
   head -1 clean00.hwe > unaff.hwe
   grep 'UNAFF' clean00.hwe >> unaff.hwe
  """
}

process generateHwePlot {
  input:
    file 'unaff.hwe' from unaff
  publishDir params.publish, overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    file 'hwe_plot.pdf'

  script:
    plotscript = script("hwe_plot_qcplink.R")
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




