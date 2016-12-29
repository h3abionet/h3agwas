#!/usr/bin/env nextflow

/*
 * Authors       :
 *
 *
 *      Shaun Aron
 *      Rob Clucas
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

import java.nio.file.Paths;
import sun.nio.fs.UnixPath;
import java.security.MessageDigest;


def helps = [ 'help' : 'help' ]
params.help = false


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

nextflowversion =getres("nextflow -v")


if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"

report = new LinkedHashMap()
repnames = ["dups","cleaned","misshet","mafpdf","snpmiss","indmisspdf","failedsex","misshetremf","diffmissP","diffmiss","pca","hwepdf","related","inpmd5","outmd5"]


repnames.each { report[it] = Channel.create() }

repmd5   = report["inpmd5"]
orepmd5 = report["outmd5"]

max_plink_cores = params.max_plink_cores 
plink_mem_req = params.plink_process_memory
other_mem_req = params.other_process_memory
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


// Checks if the file exists
checker = { fn ->
   if (fn.exists())
       return fn;
    else
       error("\n\n-----------------\nFile $fn does not exist\n\n---\n")
}



// Helper routines for phasing

/* phase takes a closure which returns the key that is used for phasing
 *  gBase returns
 *       if the parameter is a singleton
 *           the baseName if it's a file otherwise the thing itself
 *       if the paramter is a list we recursively call gBase on the first element
 */
public String gBase(fn) {
   ft = file("tmp").getClass()
   lt = [].getClass()
   if (fn.getClass() == lt) {
      res = gBase(fn[0])
   } else {
      if (fn.getClass() == ft)
      	 res = fn.baseName;
      else
         res = fn
   }
   res=res.replaceAll(/-nd.*/,"")
   return res
}	



/* combineElts is used to flatten two phased lists intelligently
 *    sometimes when we phase we have a tag value as the zeroth element
 *    of a tuple to act as a key. In this case we don't want to duplicate the key
 *    combining u and v depends on whether they are singleton values or not
 *    and if lists what the type of the zeroth element of the list is
 *    the general rule is that we remove the zeroth element of v if v is list 
 *    and the zeroth element of v is a singleton value
 */

def combineElts = { u, v ->
  ut=u.getClass()
  vt=v.getClass()
  ft = file("tmp").getClass()
  lt = [].getClass()
  st = "tmp".getClass()
  inty = 1.getClass()
  // if the two arguments are singletons return the list
  if  ( (ut in [st,inty,ft]) && vt in [st,inty,ft]) {
    return [u,v]
  }
  else
    if (ut != lt)   u = [u]
    if (vt == lt) 
       if (v[0].getClass() == ft) 
         return u+v;
       else 
         // The zeroth element of v is a key and we don't want to repeat
         return u+v[1..-1]
    else 
    return u+[v];
 
}






// Take a list of channels and phase them
def phaseAllLists = { channels ->
     start = channels[0]
     for (c in channels[1..-1]){
         start = start.phase(c) { gBase(it) } .map(combineElts)
     }
     return start;
}

// Take a list of channels and phase them
def phaseAllMap = { chanmap ->
     start = chanmap[repnames[0]]
     for (c in repnames[1..-1]){
         start = start.phase(chanmap[c]) { gBase(it) } .map(combineElts)
     }
     return start;
}




raw_ch = Channel.create()
bim_ch = Channel.create()
inpmd5ch = Channel.create()
configfile = Channel.create()

/* Get the input files -- could be a glob
 * We match the bed, bim, fam file -- order determined lexicographically
 * not by order given, we check that they exist and then 
 * send the all the files to raw_ch and just the bim file to bim_ch */
inpat = "${params.input_dir}/${params.input_pat}"


Channel
.fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true){ file -> file.baseName }  \
   .ifEmpty { error "No matching plink files" }        \
   .map { a -> [checker(a[1]), checker(a[2]), checker(a[3])] }\
   .separate(raw_ch, bim_ch, inpmd5ch,configfile) { a -> [a,a[1],a,"${workflow.projectDir}/nextflow.config"] }


// Generate MD5 sums of output files
process inMD5 {
  input:
     file plink from inpmd5ch
  output:
     file(out) into report["inpmd5"]
  echo true
  script:
       bed = plink[0]
       bim = plink[1]
       fam = plink[2]
       out  = "${plink[0].baseName}.md5"
       template "md5.py"
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
    file(inpfname) from bim_ch
  output:
    file("${base}.dups") into duplicates_ch
  script:
     print inpfname
     base = inpfname.baseName
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
    file(plinks)	\
      from raw_ch.phase(duplicates_ch){gBase(it) }.map {it.flatten()}

  publishDir params.output_dir, pattern: "$logfile", \
             overwrite:true, mode:'copy'

  output:
    set  file("${nodup}.bed"),file("${nodup}.bim"),file("${nodup}.fam")\
         into (sex_check_ch,missing_ch,het_ch,ibd_prune_ch,remove_inds_ch)
    set file("${base}.orig"), file(dups) into report["dups"]
  script:
   base    = plinks[0].baseName
   dups    = plinks[3]
   nodup   = "${base}-nd"
   logfile = "${base}-0002-dups.log"
   """
    plink --bfile $base $sexinfo --exclude $dups --make-bed --out $nodup
    wc -l ${base}.bim > ${base}.orig
    wc -l ${base}.fam >> ${base}.orig
    mv ${nodup}.log $logfile
   """
}





/* Process to identify individual discordant sex information.
 * results are put in the output directory
 */
process identifyIndivDiscSexinfo {
  memory other_mem_req
  input:
     file(plinks) from sex_check_ch

  publishDir params.output_dir, overwrite:true, mode:'copy'

  output:
     set val(base), file(logfile) into failed_sex_check
     file("${base}*.sexcheck") into report["failedsex"]
  script:
    base = plinks[0].baseName
    logfile= "${base}-failed.sexcheck"
    if (params.sexinfo_available == "true")
    """
       plink --bfile $base --check-sex  --out $base
       if grep -Rn 'PROBLEM' ${base}.sexcheck > $logfile; then
         echo 'Discordant sex info found'
       else
         echo 'No discordant sex info found'
       fi

    """
    else
     "echo 'No sex information available to check'  > $logfile"

}


// Find missingness statistics for the plink file
process calculateSampleMissing {
  memory plink_mem_req
  input:
     file(plinks) from missing_ch

  publishDir params.output_dir, overwrite:true, mode:'copy'

  output:
     file("${imiss}.imiss") into\
         (plot1_ch_miss,missing2_ch,miss_het_ch)
  script:
     base = plinks[0].baseName
     imiss= "${base}"
     """
       plink --bfile $base $sexinfo --geno 0.1 --missing --out $imiss
     """
}



// But for the  moment let's deal with heterozygosity

process calculateSampleHetrozygosity {
   memory plink_mem_req
   input:
      file(nodups) from het_ch

   publishDir params.output_dir, overwrite:true, mode:'copy'

   output:
      file("${hetf}.het") into (hetero_check_ch, plot1_ch_het)
   script:
      base = nodups[0].baseName
      hetf = "${base}"
   """
     plink --bfile $base $sexinfo --geno 0.1 --make-bed --out temp
     plink --bfile temp  $sexinfo --het  --out $hetf
   """
}



process generateMissHetPlot {
  memory other_mem_req
  errorStrategy 'ignore'

  input:
    set file(imiss), file(het) from \
        plot1_ch_miss.phase(plot1_ch_het) { gBase(it) }
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  output:
    file(pairs) into report["misshet"]

  script:
    base = imiss.baseName
    pairs   = "${base}-imiss-vs-het.pdf"
    meanhet = "${base}-meanhet_plot.pdf"
    template "miss_het_plot_qcplink.R"
}



// Find those who have too high missingness, or bad heterozygosity
process getBadIndivs_Missing_Het {
  errorStrategy 'ignore'
  memory other_mem_req
  input:
   set file(imiss), file(het) from \
       miss_het_ch.phase(hetero_check_ch) { gBase(it) }
  output:
    set val(base), file(outfname) into failed_miss_het
    file(outfname) into report["misshetremf"]
  script:
    base = het.baseName
    outfname = "${base}-fail_miss_het.txt"
    template "select_miss_het_qcplink.pl"

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
  cpus max_plink_cores
  memory plink_mem_req
  input:
    file plinks from ldreg_ch.cross(ibd_prune_ch) {true} .map { it.flatten() }
  output:
    file "${outf}.genome" into sort_ibd_ch
  script:
    base   = plinks[1].baseName
    ldreg  = plinks[0]
    outf   =  base
    if (params.high_ld_regions_fname != "")
      range = "--range --exclude $ldreg"
    else
      range =""
    """
      plink --bfile $base --threads $max_plink_cores --autosome $sexinfo $range --indep-pairwise 50 5 0.2 --out ibd
      plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --out ibd_prune
      plink --bfile $base --threads $max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --min $pi_hat --out $outf
      echo DONE
     """

}



// run script to find related individuals
//  Future - perhaps replaced with Primus
//  Currently we remove one element from each pair -- choose
//  the one with the greater missingness
process findRelatedIndiv {
  errorStrategy 'ignore'
  memory other_mem_req
  input:
     set file (missing), file (ibd_genome) from \
         missing2_ch.phase(sort_ibd_ch) { gBase(it) }
  output:
     set val(base), file(outfname) into related_indivs
     file(outfname) into report["related"]
  script:
     base = missing.baseName
     outfname = "${base}-fail_IBD.txt"
     template "run_IBD_QC_qcplink.pl"
}



process removeQCIndivs {
  memory plink_mem_req
  input:
    set val(label), file(f_miss_het), file (rel_indivs), file (f_sex_check_f),\
        file(bed), file(bim), file(fam) from\
       phaseAllLists([failed_miss_het,related_indivs,failed_sex_check,remove_inds_ch])
  script:
  output:
     file("${out}.{bed,bim,fam}") into\
         (clean00_ch1,clean00_ch2,clean00_ch3, clean00_ch4)
  script:
   base = bed.baseName
   out  = "${base}-c"
    """
     cat $f_sex_check_f $rel_indivs $f_miss_het | sort -k1 | uniq > failed_inds
     plink --bfile $base $sexinfo --remove failed_inds --make-bed --out $out
     mv failed_inds ${out}.irem
  """
}




process calculateMaf {
  memory plink_mem_req
  input:
    file(plinks) from clean00_ch1

  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.frq"

  output:
    file "${base}.frq" into maf_plot_ch

  script:
    base = plinks[0].baseName
    """
      plink --bfile $base $sexinfo  --freq --out $base
    """
}


process generateMafPlot {
  memory other_mem_req
  input:
    file frqfile from maf_plot_ch
  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  output:
    file(ofname) into report["mafpdf"]

  script:
    base    = frqfile.baseName
    ofname  = "${base}-maf_plot.pdf"
    template "maf_plot_qcplink.R"
}


// Repeat computation of missingness on QCd data
process calculateSnpMissingness {
  memory plink_mem_req
  input:
   file (plinks)  from clean00_ch2

  output:
   file ("${base}.lmiss") into clean_lmiss_plot_ch
   file ("${base}.imiss") into clean_imiss_plot_ch
  script:
   base=plinks[0].baseName
   """
     plink --bfile $base $sexinfo --missing --out $base
   """
}


process generateSnpMissingnessPlot {
  memory other_mem_req
  input:
      file(lmissf) from clean_lmiss_plot_ch

  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  output:
     file(output) into report['snpmiss']

  echo true
  script:
    input  = lmissf
    base   = lmissf.baseName
    output = "${base}-snpmiss_plot.pdf"
    template "snpmiss_plot_qcplink.R"
}


process generateIndivMissingnessPlot {
  memory other_mem_req
  input:
      file(imissf) from clean_imiss_plot_ch

  publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"

  output:
    file(output) into report["indmisspdf"]

  script:
    input  = imissf
    base   = imissf.baseName
    output = "${base}-indmiss_plot.pdf"
    template "imiss_splot.R"
}



// Find differential missingness between cases and controls; also compute HWE scores
process calculateSnpSkewStatus {
  memory plink_mem_req
  cpus max_plink_cores
  input:
    file(plinks)  from clean00_ch3
  output:
    file "${base}.missing" into clean_diff_miss_plot_ch1
    file "${base}.missing.mperm" into clean_diff_miss_ch2
    file "${base}.hwe" into hwe_scores_ch
  script:
   base = plinks[0].baseName
   """
    plink --threads ${max_plink_cores} --bfile $base $sexinfo --test-missing mperm=20000 --hardy --out $base
   """
}


process generateDifferentialMissingnessPlot {
   memory other_mem_req
   input:
     file clean_missing from clean_diff_miss_plot_ch1
   publishDir params.output_dir, overwrite:true, mode:'copy', pattern: "*.pdf"
   output:
      file output into report["diffmissP"]
   script:
       input = clean_missing
       base  = clean_missing.baseName
       output= "${base}-diff-snpmiss_plot.pdf"
       template "diffmiss_splot_qcplink.R"

 }


// Find those SNPs that have diff missingness in cases & controls
process findSnpExtremeDifferentialMissingness {
  memory other_mem_req
  input:
    file clean_missing from clean_diff_miss_ch2
  echo true
  output:
     set val(base), file(failed) into bad_snps_ch
     file(failed) into report["diffmiss"]
  script:
    cut_diff_miss=params.cut_diff_miss
    missing = clean_missing
    base     = missing.baseName
    probcol = '3'  // 4 if using the non mperm
    failed   = "${base}-failed_diffmiss.snps"
    template "select_diffmiss_qcplink.pl"
}

// Find HWE scores of each SNP
process findHWEofSNPs {
  memory other_mem_req
  input:
     file hwe from hwe_scores_ch
  output:
     file output  into unaff_hwe

  script:
    base   = hwe.baseName
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
    file output into report["hwepdf"]

  script:
    input  = unaff
    base   = unaff.baseName
    output = "${base}-hwe_plot.pdf"
    template "hwe_plot_qcplink.R"
}



process removeQCPhase1 {
  memory plink_mem_req
  input:
    file inputs  from \
         clean00_ch4\
           .phase(bad_snps_ch) { gBase(it) }\
           .map (combineElts)
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file("${output}*.{bed,bim,fam,irem,log}") into report["cleaned"]
    file("${output}*") into pca_ch
    file("${output}*.{bed,bim,fam}") into outmd5_ch
  script:
     base=inputs[0].baseName
     bad =inputs[3]
     output = "${base}-clean"
     """
     # remove really realy bad SNPs and really bad individuals
     touch temp1.irem
     plink --bfile $base $sexinfo --exclude $bad --mind 0.2 --make-bed --out temp1
     # remove bad SNPs
     plink --bfile temp1 $sexinfo --geno 0.2 --make-bed --out temp2
     # Now do final QC
     plink --bfile temp2  $sexinfo \
         --autosome \
         --maf $params.cut_maf --mind $params.cut_mind\
          --geno $params.cut_geno --hwe $params.cut_hwe \
         --make-bed --out $output 
     touch ${output}.irem
     cat temp1.irem >> ${output}.irem
  """
}

// Generate MD5 sums of output files
process outMD5 {
  input:
     file plink from outmd5_ch
  output:
     file(out) into report["outmd5"]
  echo true
  script:
       bed = plink[0]
       bim = plink[1]
       fam = plink[2]
       out  = "${plink[0].baseName}.md5"
       template "md5.py"
}




process compPCA {
   cpus max_plink_cores
   input:
      file plinks from pca_ch
   output:
    set file ("${prune}.eigenval"), file("${prune}.eigenvec") into pcares
   script:
      base = plinks[0].baseName
      prune= "${base}-prune"
     """
     plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
     plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
     plink --bfile ${prune} --pca --out $prune 
     """
}

process drawPCA {
    input:
      set file(EIGVALS), file(EIGVECS) from pcares
    output:
      file (OUTPUT) into report["pca"]
      file "rversion" into rversion
    publishDir params.output_dir, overwrite:true, mode:'copy',pattern: "*.pdf"
    script:
      base=EIGVALS.baseName
      OUTPUT="${base}-pca.pdf"
      template "drawPCA.R"

}

process produceReports {
  input:
     file results from  phaseAllMap(report)
     file rversion
     file configfile
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    file("${base}.pdf")
   script:
     println results
     orig = results[0]
     dupf = results[1]
     base = orig.baseName
     cbed = results[2]
     cbim = results[3]
     cfam = results[4]
     irem = results[5]
     ilog = results[6]
     missingvhetpdf = results[7]
     mafpdf   = results[8]
     snpmisspdf = results[9]
     indmisspdf = results[10]
     fsex        = results[11]
     misshetremf  = results[12]
     diffmisspdf  = results[13]
     diffmiss     = results[14]
     pcapdf       = results[15]
     hwepdf       = results[16]
     relf         = results[17]
     inpmd5    = results[18]
     outmd5  = results[19]
     nextflowconfig= configfile
     template "qcreport.py"
}

