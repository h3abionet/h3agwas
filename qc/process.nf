#!/usr/bin/env nextflow
include {strmem}  from '../modules/fct_groovy.nf'

/*
 * Authors       :
 *
 *      Jean-Tristan Brandenburg
 *      Scott Hazelhurst
 *
 *  On behalf of the H3ABionet Consortium and SBIMB
 *  2015-2024
 *
 *
 * Description  : 
 *
 *(C) University of the Witwatersrand, Johannesburg, 2016-2024
 *    on behalf of the H3ABioNet Consortium
 *This is licensed under the MIT Licence. See the "LICENSE" file for details
 */

//---- General definitions --------------------------------------------------//

                                                                                
                                                                                
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
              




def checkSampleSheet(fname)  {
  if (workflow.profile == "awsbatch") return;
  if (fname.contains("s3://") )return;
  if (fname.contains("az://") ) return;
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




//This method first checks that the data file has the stated column 
// If so, it creates a channel for it
// NB: if the file is in S3 we cannot do the test since Groovy does not
// allow us to access the file directly
def fileheader_create_ch = { parm, parm_name, col_name ->
  if (parm.toString().contains("s3://")) {
    println "The file <$parm> is in S3 so we cannot do a pre-check";
    return Channel.fromPath(parm, check);
  }
  if (parm.toString().contains("az://")) {
    println "The file <$parm> is in Azure so we cannot do a pre-check";
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




//---- Modification of variables for pipeline -------------------------------//


// Generate MD5 sums of output files
process MD5 {
  input:
     path(plink) 
  output:
     path(out) 
  script:
       bed = plink[0]
       bim = plink[1]
       fam = plink[2]
       out  = "${plink[0].baseName}.md5"
       template "md5.py"
}



  process sampleSheet {
    input:
       tuple path(sheet), val(balise), val(idpat)
    output:
     path("poorgc10.lst"), emit: poorsgc10
     path("plates") , emit: plates
    script:
     println balise
     if(balise.toString()=="1"){
     """
       mkdir -p plates
       sampleqc.py $sheet ${params.gc10} "${idpat}"  poorgc10.lst plates/crgc10.tex
      """
     }else {
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
  memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  publishDir "${params.output_dir}/snps/duplicate_marker", pattern: "*dups", \
             overwrite:true, mode:'copy'
  input:
    path(inpfname) 
    val(remove_on_bp)
  output:
    path("${base}.dups") 
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
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
   tuple path(bed), path(bim), path(fam) 
   path(dups) 
   val(extrasexinfo)
   val(sexinfo)
  output:
    tuple path("${nodup}.bed"),path("${nodup}.bim"),path("${nodup}.fam"), emit : plink
    tuple path("${base}.orig"), path(dups), emit : dups 
    path ("${nodup}.lmiss"), emit : lmiss
    path ("${nodup}.imiss"), emit : imiss
  script:
   base    = bed.baseName
   nodup   = "${base}-nd"
   K = "--keep-allele-order"
   """
    plink $K --bfile $base $sexinfo $extrasexinfo --exclude $dups --missing --make-bed --out $nodup
    wc -l ${base}.bim > ${base}.orig
    wc -l ${base}.fam >> ${base}.orig
   """
}


/*process to check if X */
process countChr{
 input :
  tuple path(bed), path(bim), path(fam)
  val(chro)
 output :
   stdout
 script :
 """
 awk '{if(\$1==${chro})print \$1}' $bim |wc -l
 """
}




process clean_x {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }

   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input :   
    tuple path(bed), path(bim), path(fam)  
    val(cxy) 
    val(cxx) 
  output :
    tuple path("${baseout}.bed"), path("${baseout}.bim"), path("${baseout}.fam") 
  script :
   base=bed.baseName
   cxy=cxy.toInteger()
   cxx=cxx.toInteger()
    if(cxy==0 && cxx>0){
        baseout=base+'_splx'
        """
        plink -bfile $base   --make-bed --out $baseout --split-x $params.build_genome 'no-fail' --keep-allele-order
        """
    }else{
       baseout=base
       """
       echo $baseout
       """
    }
}

process getX {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

     input:
       path(plink) 
      output:
       path("X*") 
      script:
      base = plink[0].baseName
      """
	if [[ `grep  "^${params.chrxx_plink}" *bim`  ]];  then
	   plink --bfile $base --chr ${params.chrxx_plink} --geno 0.04 --make-bed --out X
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


   // input : bim contained x
process analyseX {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3
     input:
       path(xchr) 
     output:
       path(out) 
     script:
        f_hi_female = params.f_hi_female
        f_lo_male = params.f_lo_male
        missingness = [0.01,0.03,0.05]  // this is used by one of the templates
	x = xchr[0].baseName
	out = "x.pkl"
	template "xCheck.py"
}

   

/* Process to identify individual discordant sex information.
 * results are put in the output directory
 * Also does HWE
 */
process identifyIndivDiscSexinfo {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
     path(plinks) 

  publishDir "${params.output_dir}/samples/sexinfo", overwrite:true, mode:'copy'
  output:
     path(logfile), emit : log
     tuple path(imiss), path(lmiss),path(sexcheck_report), emit : stat
     path("${base}.hwe"), emit : hwe 
  errorStrategy { task.exitStatus in [0,1] ? 'ignore' : 'terminate' }
  script:
    base = plinks[0].baseName
    logfile= "${base}.badsex"
    sexcheck_report = "${base}.sexcheck"
    imiss  = "${base}.imiss"
    lmiss  = "${base}.lmiss"
    K = "--keep-allele-order"
    if (params.sexinfo_available == true)
    """
       plink $K --bfile $base --hardy --check-sex $params.f_hi_female $params.f_lo_male --missing  --out $base
       head -n 1 ${base}.sexcheck > $logfile
       grep  'PROBLEM' ${base}.sexcheck >> $logfile
    """
    else
     """
     plink --bfile $base  --hardy --missing  --out $base $K
     echo 'FID IID STATUS' > $sexcheck_report
     echo 'No sex information'  > $logfile
     """
}



process generateSnpMissingnessPlot {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
      path(lmissf)
  publishDir "${params.output_dir}/snps/missingness", overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
     path(output) 
  script:
    input  = lmissf
    base   = lmissf.baseName
    label  = "SNPs"
    output = "${base}-snpmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}


process generateIndivMissingnessPlot {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
      path(imissf) 
  publishDir "${params.output_dir}/samples/missingness", overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    path(output) 
  script:
    input  = imissf
    base   = imissf.baseName
    label  = "samples"
    output = "${base}-indmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}
 
process getInitMAF {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
     path(plink) 
  output:
     path("${newbase}.frq") 
  script:
    base = plink[0].baseName
    newbase = base.replace(".","_")
    """
    plink --bfile $base --freq --out $newbase
    """
}


process showInitMAF {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
     path(freq) 
  output:
     tuple path("${base}.pdf"), path("${base}.tex") 
  script:
    base = freq.baseName+"-initmaf"
    base = base.replace(".","_")
    template "showmaf.py"
}

process showHWEStats {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
     path(hwe) 
  output:
     tuple path("${base}.pdf"), path("${base}-qq.pdf"), path("${base}.tex") 
  script:
    base = hwe.baseName+"-inithwe"
    base = base.replace(".","_")
    template "showhwe.py"
}


process removeQCPhase1 {
  memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    tuple path(bed), path(bim), path(fam) 
    val(sexinfo)
  publishDir "${params.output_dir}/phase1/", overwrite:true, mode:'copy'
  output:
    path("${output}*.{bed,bim,fam}"), emit : plink
    tuple path("qc1.out"), path("${output}.irem"), emit : log
  script:
     base=bed.baseName
     output = "${base}-c".replace(".","_")
     K = "--keep-allele-order"
     """
     # remove really realy bad SNPs and really bad individuals
     plink $K ${params.autosome_plink} --bfile $base $sexinfo --mind 0.1 --geno 0.1 --make-bed --out temp1
     plink $K --bfile temp1  $sexinfo --mind $params.cut_mind --make-bed --out temp2
     /bin/rm temp1.{bim,fam,bed}
     plink $K --bfile temp2  $sexinfo --geno $params.cut_geno --make-bed --out temp3
     /bin/rm temp2.{bim,fam,bed}
     plink $K --bfile temp3  $sexinfo --maf $params.cut_maf --make-bed --out temp4
     /bin/rm temp3.{bim,fam,bed}
     plink $K --bfile temp4  $sexinfo --hwe $params.cut_hwe --make-bed  --out $output
     /bin/rm temp4.{bim,fam,bed}
     cat *log > logfile
     touch tmp.irem
     cat *.irem > ${output}.irem
     qc1logextract.py logfile ${output}.irem > qc1.out     
  """
}


// We do PCA on qc2 data because relatively few SNPs and individuals will be removed later and
// this is an expensive operation so we start early. Similarly for computing relatedness
process compPCA {
  cpus params.max_plink_cores
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
   input:
      path(plinks) 
   output:
      tuple path("${prune}.eigenval"), path("${prune}.eigenvec"), emit:eigen
      tuple path ("${prune}.bed"), path("${prune}.bim"), path("${prune}.fam"), emit : plink_prune
   publishDir "${params.output_dir}/pca", overwrite:true, mode:'copy',pattern: "${prune}*"
   script:
      base = plinks[0].baseName
      prune= "${base}-prune".replace(".","_")
     """
     plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
     plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
     /bin/rm check*
     plink --bfile ${prune} --pca --out $prune 
     """
}

process drawPCA {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

    input:
      tuple path(eigvals), path(eigvecs) 
      path(cc) 
      val(col) 
      val(diffpheno) 

    output:
      tuple  path ("eigenvalue.pdf"), path(output) 
    publishDir "${params.output_dir}/pca/", overwrite:true, mode:'copy',pattern: "*.pdf"
    script:
      base=eigvals.baseName
      cc_fname = params.case_control
      // also relies on "col" defined above
      output="${base}-pca".replace(".","_")+".pdf"
      template "drawPCA.py"

}




  process pruneForIBDLD {
    cpus params.max_plink_cores
    memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 3

    input:
      path(plinks)
      path(ldreg)  
      val(sexinfo)
      val(pi_hat)
    publishDir "${params.output_dir}/samples/ibd", overwrite:true, mode:'copy'
    output:
      path("${outf}.genome") 
    script:
      base   = plinks[0].baseName
      outf   =  base.replace(".","_")
      range = ldreg.name != 'NO_FILE' ?  " --exclude range $ldreg" : ""
      """
       plink --bfile $base --threads $params.max_plink_cores --autosome $sexinfo $range --indep-pairwise 60 5 0.2 --out ibd
       plink --bfile $base --threads $params.max_plink_cores --autosome $sexinfo --extract ibd.prune.in --genome --min $pi_hat --out $outf
       echo LD
       """
  }

// run script to find a tuple of individuals we can remove to ensure no relatedness
//  Future - perhaps replaced with Primus
process findRelatedIndiv {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
     path(missing) 
     path(ibd_genome) 
  output:
     path(outfname) 
  publishDir "${params.output_dir}/samples/relatdness", overwrite:true, mode:'copy'
  script:
     super_pi_hat = params.super_pi_hat
     pi_hat = params.pi_hat
     base = missing.baseName
     outfname = "${base}-fail_IBD".replace(".","_")+".txt"
     template "removeRelInds.py"
}



process calculateSampleHeterozygosity {
   memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

   input:
      path(nodups) 
      val(sexinfo)

   publishDir "${params.output_dir}/samples/heterozygoty", overwrite:true, mode:'copy'
   output:
      tuple path("${hetf}.het"), path("${hetf}.imiss"), emit: hetmiss
      path("${hetf}.imiss"), emit: imiss
   script:
      base = nodups[0].baseName
      hetf = "${base}".replace(".","_")
   """
     plink --bfile $base  $sexinfo --het --missing  --out $hetf
   """
}



process generateMissHetPlot {
   memory { strmem(params.low_memory) + 5.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3
  input:
    tuple path(het), path(imiss) 
  publishDir "${params.output_dir}/samples/heterozygoty", overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    path(output)
  script:
    base = imiss.baseName
    output  = "${base}-imiss-vs-het".replace(".","_")+".pdf"
    template "missHetPlot.py"
}



// Find those who have bad heterozygosity
process getBadIndivsMissingHet {
   memory { strmem(params.low_memory) + 5.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
    tuple path(het), path(imiss) 
  output:
    path(outfname) 
  publishDir "${params.output_dir}/samples/heterozygoty", overwrite:true, mode:'copy', pattern: "*.txt"
  script:
    base = het.baseName
    outfname = "${base}-fail_het".replace(".","_")+".txt"
    template "select_miss_het_qcplink.py"
}




process removeQCIndivs {
   memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
    path(f_miss_het)     
    path(rel_indivs)     
    path (f_sex_check_f) 
    path (poorgc)       
    tuple path(bed), path(bim), path(fam) 
    val(sexinfo)
  output:
     path("${out}.{bed,bim,fam}"), emit: plink
     path("${out}.fam"), emit: fam
  script:
   base = bed.baseName
   out  = "${base}-c".replace(".","_")
    K = "--keep-allele-order"
    """
     cat $f_sex_check_f $rel_indivs $f_miss_het $poorgc | sort -k1 | uniq > failed_inds
     plink $K --bfile $base $sexinfo --remove failed_inds --make-bed --out $out
     mv failed_inds ${out}.irem
  """
}


mperm_header=" CHR                               SNP         EMP1         EMP2 "

// Find differential missingness between cases and controls; also compute HWE scores
process calculateSnpSkewStatus {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  cpus params.max_plink_cores
  input:
   path(plinks) 
   path(phe)
   val(sexinfo)
   val(diffpheno)
  output:
    path("${base}.missing"), emit: missing 
    path(mperm) , emit: perm
    path("${base}.hwe") , emit :hwe
  script:
   base  = plinks[0].baseName
   out   = base.replace(".","_")
   mperm = "${base}.missing.mperm"
   """
    cp $phe cc.phe
    plink --threads ${params.max_plink_cores} ${params.autosome_plink} --bfile $base $sexinfo $diffpheno --test-missing mperm=10000 --hardy --out $out
    if ! [ -e $mperm ]; then
       echo "$mperm_header" > $mperm
    fi

   """
}


process generateDifferentialMissingnessPlot {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

   input:
     path(clean_missing) 
   publishDir "${params.output_dir}/snps/missingness", overwrite:true, mode:'copy', pattern: "*.pdf"
   output:
      path(output) 
   script:
       input = clean_missing
       base  = clean_missing.baseName.replace(".","_").replace("-nd","")
       output= "${base}-diff-snpmiss_plot.pdf"
       template "diffMiss.py"

 }


// Find those SNPs that have diff missingness in cases & controls
process findSnpExtremeDifferentialMissingness {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3
  input:
    path(clean_missing) 
  echo true
  output:
     tuple val(base), path(failed), emit : base_failed
     path(failed), emit : failed
  script:
    cut_diff_miss=params.cut_diff_miss
    missing = clean_missing
    base     = missing.baseName.replace("-.*","").replace(".","_")
    probcol = 'EMP2'  // need to change if we don't use mperm
    failed   = "${base}-failed_diffmiss.snps"
    template "select_diffmiss_qcplink.py"
}

process removeSkewSnps {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    path(plinks) 
    path(failed) 
    val(sexinfo)
  //publishDir "${params.output_dir}/", overwrite:true, mode:'copy'  
  output:
    tuple path("${output}.bed"), path("${output}.bim"), path("${output}.fam"), path("${output}.log") 
  script:
  base = plinks[0].baseName
  output = params.output.replace(".","_")
  K = "--keep-allele-order"
  """
  plink $K --bfile $base $sexinfo --exclude $failed --make-bed --out $output
  """
}


process splitX {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input :
   tuple path(bed), path(bim), path(fam) 
   path(listind) 
  output :
   tuple  path("${outputx}.bed"),path("${outputx}.bim"),path("${outputx}.fam") 
  script :
   base = bed.baseName
   outputx="{params.output}_x"
   """
   plink -bfile $base --chr ${params.chrxx_plink} --make-bed --out $outputx --keep-allele-order --keep $listind
   """
 }

 process cleanPlink_x {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
  input :
   tuple path(bed), path(bim), path(fam) 
  output :
   tuple  path("${outputfem}.bed"),path("${outputfem}.bim"),path("${outputfem}.fam"), emit : female_x
   tuple  path("${outputmale}.bed"),path("${outputmale}.bim"),path("${outputmale}.fam"), emit :male_x
   tuple path("list_female"), path("list_male"), emit :listind
  script  :
   base = bed.baseName
   outputx="${params.output}_x"
   outputfem="${params.output}_x_female"
   outputmale="${params.output}_x_male"
   """
   plink -bfile $base --chr ${params.chrxx_plink} -make-bed --out $outputx
   awk '{if(\$5==2)print \$1"\\t"\$2}' $outputx".fam"  > list_female
   plink -bfile $base --keep list_female --chr ${params.chrxx_plink} -make-bed --out $outputfem
   awk '{if(\$5==1)print \$1"\\t"\$2}' $outputx".fam"  > list_male
   plink -bfile $base --keep list_male --chr ${params.chrxx_plink} -make-bed --out $outputmale --set-hh-missing
   """ 
 }

process computed_stat_female_x {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input :
   tuple path(bed), path(bim), path(fam) 
  publishDir "${params.output_dir}/qcX/female", overwrite:true, mode:'copy'
  output :
    tuple path("${output}.frq"), path("${output}.hwe"), path("${output}.imiss"),  path("${output}.lmiss") 
  script :
    base = bed.baseName
    output=params.output+"_femalestat"
    """
    plink -bfile $base --missing  --freq  --hardy  -out $output
    """
  }
 
process computed_stat_male_x {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input :
   tuple path(bed), path(bim), path(fam) 
  publishDir "${params.output_dir}/qcX/male", overwrite:true, mode:'copy'
  output :
    tuple path("${output}.frq"), path("${output}.imiss"),  path("${output}.lmiss")
  script :
    base = bed.baseName
    output=params.output+"_malestat"
    """
    plink -bfile $base --missing  --freq    -out $output
    """
  }

process report_export_x {
    input :
      tuple path(frqmale),  path(maleimiss),  path(malelmiss) 
      tuple path(frqfem), path(femhwe), path(femimiss),  path(femlmiss) 
      publishDir "${params.output_dir}/qcX/"
      output :
	     path("${output}.tex"), emit : report
	     path("${output}.in"), emit : listsnps
	     path("${output}_resume.csv") , emit : resume
      script :
	      basemal=frqmale.baseName
	      basefem=frqfem.baseName
	      output=params.output+"_xfilter"
	      """
	      stats_x.py --base_female  $basefem --base_male $basemal --out $output --maf_female ${params.cut_maf_xfemale} --maf_male ${params.cut_maf_xmale}  --miss_male ${params.cut_miss_xmale} --miss_female ${params.cut_miss_xfemale} --diff_miss ${params.cut_diffmiss_x}
	      """
}

 process cleanandmerge_x {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
 input :
   tuple path(bed), path(bim), path(fam), path(log) 
	   tuple path(bedx), path(bimx), path(famx) 
	   path(snpsx) 
	  publishDir "${params.output_dir}/", overwrite:true, mode:'copy'
  output :
	   tuple  path("${output}.bed"),path("${output}.bim"),path("${output}.fam"), path("${output}.log") 

  script :  
	    basex=bedx.baseName
	    base=bed.baseName
	    output=params.output+"_withx"
	    """
	    plink -bfile $basex --extract $snpsx --keep-allele-order  -out cleanx --make-bed
	    plink -bfile $base --bmerge cleanx  --keep-allele-order  -out $output --make-bed
	    """
 }

  process report_export_x_tmp {
	   output :
	     path("noreport") 
	   script :
	      """
	      echo "no x analyse " > noreport
	      """
}


process clean_y {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
  input :
	   tuple path(bed), path(bim), path(fam)
           path(listind) 
           val(county)
  publishDir "${params.output_dir}/qcY/", overwrite:true, mode:'copy'
  output :
            tuple path("${outputy}.bed"), path("${outputy}.bim"), path("${outputy}.fam"), emit : plky_qc
            tuple path("${outputy}.lmiss"),  path("${outputy}.imiss"), path("${outputy}.frq"), emit : stat_qc
  script :
      base=bed.baseName 
      outputy=params.output+"_y"
      county = county.toInteger()
      if(county > 0)
           """
           awk '{if(\$5==1)print \$1"\\t"\$2}' $base".fam"  > list_male
           plink -bfile $base --chr ${params.chry_plink} -make-bed --out tmp1 --keep $listind #list_male --set-hh-missing  
           plink -bfile tmp1 --keep list_male --set-hh-missing  --out tmp --make-bed
           plink -bfile tmp --missing  --freq  --out $outputy
           plink -bfile tmp --maf ${params.cut_maf_y} --geno ${params.cut_miss_y} --make-bed --out tmpy
           awk '{print \$2}' tmpy".bim" > rslist
           plink -bfile $base --extract rslist --make-bed --set-hh-missing -out $outputy --keep-allele-order
           """
      else 
          """
          touch "${outputy}.bed"  "${outputy}.bim" "${outputy}.fam" "${outputy}.lmiss" "${outputy}.imiss" "${outputy}.frq"
          """
}

process cleanandmerge_y {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
          input :
           tuple path(bed), path(bim), path(fam), path(log) 
           tuple path(bedy), path(bimy), path(famy) 
           val(county)
          publishDir "${params.output_dir}/qcY/", overwrite:true, mode:'copy'
          output :
           tuple  path("${output}.bed"),path("${output}.bim"),path("${output}.fam"), path("${output}.log") 
          script :
            basey=bedy.baseName
            base=bed.baseName
            output=params.output
            county = county.toInteger()
            if(county>0)
              """
              plink -bfile $base --bmerge $basey --keep-allele-order  -out $output --make-bed
              """
            else 
              """
              cp ${base}.bed ${output}.bed
              cp ${base}.bim ${output}.bim
              cp ${base}.fam ${output}.fam
              cp ${base}.log ${output}.log
              """
}

process build_reporty{
           input :
            tuple path(lmiss),  path(imiss), path(frq)
            val(county)
          publishDir "${params.output_dir}/qcY/", overwrite:true, mode:'copy'
          output :
             path("${output}.tex") , emit : report
             path("${output}.in"), emit : listsnps
             path("${output}_resume.csv"), emit : resume
          script :
             output=lmiss.baseName
             county = county.toInteger()
             if(county>0)
               """
               stats_y.py --base $output --out $output --maf ${params.cut_maf_y} --miss ${params.cut_miss_y}
               """
            else 
              """                                                               
              echo "No Y analyse " > $output".tex"
              touch ${output}.in
              touch ${output}_resume.csv
              """   
}







process calculateMaf {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3

  input:
    tuple  path(bed), path(bim), path(fam), path(log) 
    val(sexinfo)

  publishDir "${params.output_dir}/snps/maf", overwrite:true, mode:'copy', pattern: "*.frq"

  output:
    path("${base}.frq") 

  script:
    base = bed.baseName
    out  = base.replace(".","_")
    """
      plink --bfile $base $sexinfo  --freq --out $out
    """
}


process generateMafPlot {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
    path(input) 
  publishDir "${params.output_dir}/snps/maf", overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    path(output) 
  script:
    base    = input.baseName
    output  = "${base}-maf_plot.pdf"
    template "mafplot.py"
}




// Find HWE scores of each SNP
process findHWEofSNPs {
   memory { strmem(params.low_memory) + 1.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3
  input:
     path(hwe) 
  output:
     path(output)

  script:
    base   = hwe.baseName.replace(".","_")
    output = "${base}-unaff.hwe"
    """
      head -1 $hwe > $output
      grep 'UNAFF' $hwe >> $output
    """
}

process generateHwePlot {
  memory { params.low_memory * task.attempt }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
  maxRetries 3
  input:
    path(unaff) 
  publishDir "${params.output_dir}/snps/hwe", overwrite:true, mode:'copy', pattern: "*.pdf"
  output:
    path(output) 
  script:
    input  = unaff
    base   = unaff.baseName.replace(".","_")
    output = "${base}-hwe_plot.pdf"
    template "hweplot.py"
}




process batchProc {
  memory { strmem(params.plink_mem_req) + 5.GB * (task.attempt -1) }
  errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3

  input:
    tuple path(eigenval), path(eigenvec)
    tuple path(imiss), path(lmiss), path(sexcheck_report) 
    path(pheno), name : 'pheno.phe'
    path(batch), name: 'batch.batch'
    path(genome) 
    path(pkl)   
    path(rem_indivs) 
    val(extrasexinfo)
  publishDir "${params.output_dir}/batch", pattern: "*{csv,pdf}", \
             overwrite:true, mode:'copy'
  output:
      path("${base}-batch.tex"),      emit : report_batch_report_ch
      tuple path("*.csv"), path("*pdf"), emit :report_batch_aux_ch // need to stage
  script:
    phenotype = pheno
    batch = batch
    base = eigenval.baseName
    batch_col = params.batch_col
    f_hi_female = params.f_hi_female                                        
    f_lo_male = params.f_lo_male           
    pheno_col = params.pheno_col
    pi_hat = params.pi_hat
    super_pi_hat=params.super_pi_hat
    template "batchReport.py"
}



//repnames = ["dups","cleaned","misshet","mafpdf","snpmiss","indmisspdf","failedsex","misshetremf","diffmissP","diffmiss","pca","hwepdf","related","inpmd5","outmd5","batch"]



process produceReports {
  label 'latex'
  input:
    tuple path(orig), file (dupf) 
    tuple path(cbed), path(cbim), path(cfam),path(ilog) 
    path(missingvhetpdf) 
    path(mafpdf)         
    path(snpmisspdf)    
    path(indmisspdf) 
    path(fsex)      
    path(misshetremf)
    path(diffmisspdf)
    path(diffmiss) 
    tuple path(eigenvalpdf),path(pcapdf)
    path(hwepdf)
    path(rel_indivs) 
    path(inpmd5)
    path(outmd5)
    tuple path(initmafpdf), path(initmaftex)
    tuple path(inithwepdf), path(inithweqqpdf), path(inithwetex)
    tuple path(qc1), path(irem) 
    path(batch_tex)
    path(poorgc)  
    tuple path(bpdfs), path(bcsvs)
    path(qcx)
    path(qcy) 
  publishDir params.output_dir, overwrite:true, mode:'copy'
  output:
    path("${base}.pdf") 
   script:
     base = params.output
     config_text = getConfig()
     nextflowversion =nextflow.version                                               
     if (workflow.repository)                                                        
       wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
     else                                                                            
     wflowversion="A local copy of the workflow was used"       
     template "qcreport.py"
}

/*


final_ch.subscribe { b=it.baseName; 
  new File("emptyZ0batch.txt").delete();
  new File("emptyZ0pheno.txt").delete();
  new File("xxemptyZ0pheno.txt").delete();
  println "The output report is called ${params.output_dir}/${b}.pdf"}
*/
