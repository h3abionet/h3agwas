// Generate MD5 sums of output files
plink_mem_reqq = params.plink_mem_reqq 
other_mem_reqq = params.other_mem_reqq
max_plink_cores = params.max_plink_coress

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
params.batch = "/home/wilson/sadacc-workflows/projects/testpca/input/sample-batch-site.phe"
params.phenotype="/home/wilson/sadacc-workflows/projects/testpca/input/sampleA.phe"
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

fname = params.batch

win = params.batch

def checkColumnHeader(win, columns) {
  if (workflow.profile == "awsbatch") return;
  if (fname.toString().contains("s3://")) return;
  if (nullfile.contains(fname)) return;
  new File(fname).withReader { line = it.readLine().tokenize() }  
  problem = false; 
  columns.each { col -> 
    if (! line.contains(col) ) {
      println "The file <$win> does not contain the column <$col>";
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

// println "The batch file is ${params.batch}"

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

// params.each { parm ->
//   if (! allowed_params.contains(parm.key)) {
//     	println "Check $parm  ************** is it a valid parameter -- are you using one rather than two - signs or vice-versa";
//   }
// }

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

// // This method first checks that the data file has the stated column 
// // If so, it creates a channel for it
// // NB: if the file is in S3 we cannot do the test since Groovy does not
// // allow us to access the file directly
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
  Channel.fromPath(ccfile).set { cc_ch }
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
  cc_ch  = Channel.value.into("none").set { cc_ch }
}

phenotype_ch = getSubChannel(params.phenotype,"pheno",params.pheno_col)
batch_ch     = getSubChannel(params.batch,"batch",params.batch_col)

// //---- Modification of variables for pipeline -------------------------------//


// /* Define the command to add for plink depending on whether sexinfo is
//  * available or not.
//  */

if ( nullfile.contains(params.sexinfo_available) ) {
  sexinfo = "--allow-no-sex"
  extrasexinfo = ""
//   println "Sexinfo not available, command --allow-no-sex\n"
} else {
  sexinfo = ""
  extrasexinfo = "--must-have-sex"
//   println "Sexinfo available command"
}

/* Get the input files -- could be a glob
 * We match the bed, bim, fam file -- order determined lexicographically
 * not by order given, we check that they exist and then 
 * send the all the files to raw_ch and just the bim file to bim_ch */
inpat = "${params.input_dir}"

if (inpat.contains("s3://")) {
  print "Here"
  ob = { it -> return it}
} else {
  ob = checker
}

Channel
   .fromFilePairs("${inpat}.{bed,bim,fam}",size:3, flat : true)
    { file -> file.baseName }  \
      .ifEmpty { error "No matching plink files" }        \
      .map { a -> [a[1], a[2], a[3]] }\
      .multiMap {  it ->
         raw_ch: it
         bed_ch: it[0]
         bim_ch: it[1] 
         inpmd5ch : it
       }.set {checked_input}

bim_ch = Channel.fromPath("/home/wilson/sadacc-workflows/projects/testpca/input/sampleA.bim")

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







/* ---------------------    Pipeline    ----------------------*/







process inMD5 {
  echo true
  tag "md5"

  input:
     path plink

  output:
     path out

  script:
       bed = plink[0]
       bim = plink[1]
       fam = plink[2]
       out  = "${plink[0].baseName}.md5"
       template "md5.py"
}

process sampleSheet {
  tag "Importing samplesheets"

    input:
     path sheet

    output:
     path "poorgc10.lst"
     path "plates"

    script:
     """
       mkdir -p plates
       sampleqc.py $sheet ${params.gc10} "${idpat}"  poorgc10.lst plates/crgc10.tex
      """
}

process noSampleSheet {
    tag "Importing samplesheets"

    output:
     path "poorgc10.lst"
     path "plates"

    script:
     """
       mkdir -p plates
       sampleqc.py 0 0 0 poorgc10.lst plates/crgc10.tex
      """
}

process getDuplicateMarkers {
  tag "Getting duplicate markers"

  memory other_mem_reqq
  publishDir "${params.output_dir}/DuplicateMarkers", pattern: "*dups", \
             overwrite:true, mode:'copy'

  input:
    path inpfname

  output:
    path "${base}.dups"

  script:
     base     = inpfname.baseName
     outfname = "${base}.dups"
     template "dups.py"
}

process removeDuplicateSNPs {
  tag "Removing duplicate snps"

  memory plink_mem_reqq

  input:
   path inpfname
   path dups

  output:
    path "${nodup}.bed"
    path "${nodup}.bim"
    path "${nodup}.fam"
    path "${base}.orig"
    path dups
    path "${nodup}.lmiss"
    path "${nodup}.imiss"

  script:
    base    = inpfname[0].baseName
    nodup   = "${base}-nd"

    """ echo ${nodup}"""

    """
    plink $K --bfile $base $sexinfo $extrasexinfo --exclude $dups --missing --make-bed --out $nodup
    wc -l ${base}.bim > ${base}.orig
    wc -l ${base}.fam >> ${base}.orig
    """
}

missingness = [0.01,0.03,0.05]

process getX {
    tag "X chromosomes"

    memory other_mem_req

    input:
       path bed
       path bim
       path fam
       path orig
       path dup
       path lmiss
       path imiss

    output:
       path "X*" 

    script:

      base = bed.baseName

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
       path xchr

     output:
       path out

     script:
       x = xchr[0].baseName
       out = "x.pkl"
       template "xCheck.py"
}

process identifyIndivDiscSexinfo {

  memory plink_mem_req
  publishDir "${params.output_dir}/IndivDiscSexinfo", overwrite:true, mode:'copy'

  input:
     path bed
     path bim
     path fam
     path orig 
     path dup
     path lmiss
     path imiss

  output:
     path logfile
     path sexcheck_report
     path "${base}.hwe" 
     path imiss
     path lmiss
  // validExitStatus 0, 1
  '''
  trap 'if [[ $? == expected_error ]]; then echo OK; exit 0; fi' EXIT
  your_command_that_may_fail
  '''

  script:
    base = bed.baseName
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

process identifyIndivDiscSexinfo_2 {

  memory plink_mem_req
  publishDir "${params.output_dir}/IndivDiscSexinfo", overwrite:true, mode:'copy'

  input:
     path bed
     path bim
     path fam
     path orig 
     path dup
     path lmiss
     path imiss

  output:
     path logfile
     path sexcheck_report
     path "${base}.hwe" 
  // validExitStatus 0, 1
  '''
  trap 'if [[ $? == expected_error ]]; then echo OK; exit 0; fi' EXIT
  your_command_that_may_fail
  '''

  script:
    base = bed.baseName
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
  tag "Missingness plot"

  memory other_mem_req
  publishDir "${params.output_dir}/miss", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
     path bed
     path bim
     path fam
     path orig 
     path dup
     path lmissf
     path imiss

  output:
     path output

  script:
    input  = lmissf
    base   = lmissf.baseName
    label  = "SNPs"
    output = "${base}-snpmiss_plot".replace(".","_")+".pdf"
    template "missPlot.py"
}

process generateIndivMissingnessPlot {

  memory other_mem_req
  publishDir "${params.output_dir}/miss", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
     path bed
     path bim
     path fam
     path orig 
     path dup
     path lmissf
     path imissf
  
  output:
    path output
  
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
     path bed
     path bim
     path fam
     path orig 
     path dup
     path lmissf
     path imissf

  output:
     path "${newbase}.frq"

  script:
    base = bed.baseName
    newbase = base.replace(".","_")

    """
    plink --bfile $base --freq --out $newbase
    """
}

process showInitMAF {

  memory other_mem_req
  publishDir "${params.output_dir}/InitMAF", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
     path freq

  output:
     path "${base}.pdf"
     path "${base}.tex"

  script:
    base = freq.baseName+"-initmaf"
    base = base.replace(".","_")
    template "showmaf.py"
}

process showHWEStats {

  memory other_mem_req
  publishDir "${params.output_dir}/HWEStats", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
     path logfile
     path sexcheck_report
     path hwe
     path imiss
     path lmiss

  output:
     path "${base}.pdf"
     path "${base}-qq.pdf"
     path "${base}.tex"

  script:
    base = hwe.baseName+"-inithwe"
    base = base.replace(".","_")
    template "showhwe.py"
}

process removeQCPhase1 {
  tag "Removing QC Phase 1"

  memory plink_mem_req
  publishDir "${params.output_dir}/QCPhase1", overwrite:true, mode:'copy'

  input:
    path bed
    path bim
    path fam
    path orig 
    path dup
    path lmissf
    path imissf

  output:
    path "${output}*.bed"
    path "${output}*.bim"
    path "${output}*.fam"
    path "qc1.out"
    path "${output}.irem"
  
  script:
     base = bed.baseName
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
   publishDir "${params.output_dir}/pca", overwrite:true, mode:'copy',pattern: "${prune}*"

   input:
      path bed
      path bim
      path fam
      path qc1
      path irem

   output:
      path "${prune}.eigenval"
      path "${prune}.eigenvec"
      path "${prune}.bed"
      path "${prune}.bim"
      path "${prune}.fam"
   
   script:
      base = bed.baseName
      prune= "${base}-prune".replace(".","_")
      """
      plink --bfile ${base} --indep-pairwise 100 20 0.2 --out check
      plink --bfile ${base} --extract check.prune.in --make-bed --out $prune
      plink --bfile ${prune} --pca --out $prune 
      """
}

process drawPCA {

    memory other_mem_req
    publishDir "${params.output_dir}/draw_pca", overwrite:true, mode:'copy',pattern: "*.pdf"

    input:
      path eigvals
      path eigvecs
      path bed
      path bim
      path fam
      path cc

    output:
      path "eigenvalue.pdf"
      path output
  
    script:
      base = eigvals.baseName
      cc_fname = params.case_control

      output="${base}-pca".replace(".","_")+".pdf"
      template "drawPCA.py"
}

process pruneForIBD {

  cpus max_plink_cores
  memory plink_mem_req
  publishDir "${params.output_dir}/pruneForIBD", overwrite:true, mode:'copy'

  input:
    path bed
    path bim
    path fam
    path qc1
    path irem

  output:
    path "${outf}.genome"

  script:
    base   = bed.baseName
    outf   =  base.replace(".","_")

    if (params.high_ld_regions_fname != "") {
      ldreg_ch=Channel.fromPath(params.high_ld_regions_fname)
      range = " --exclude range $ldreg_ch"
    }
    else{
      range =""
    }

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
  publishDir "${params.output_dir}/findRelatedIndiv", overwrite:true, mode:'copy'

  input:
    path bed
    path bim
    path fam
    path orig 
    path dup
    path lmissf
    path missing
    path ibd_genome

  output:
     path outfname

  script:
     base = missing.baseName
     outfname = "${base}-fail_IBD".replace(".","_")+".txt"
     template "removeRelInds.py"
}

process calculateSampleHeterozygosity {

   memory plink_mem_req
   publishDir "${params.output_dir}/Heterozygosity", overwrite:true, mode:'copy'

   input:
      path bed
      path bim
      path fam
      path qc1
      path irem

   output:
      path "${hetf}.het"
      path "${hetf}.imiss"

   script:
      base = bed.baseName
      hetf = "${base}".replace(".","_")
      """
        plink --bfile $base  $sexinfo --het --missing  --out $hetf
      """
}

process generateMissHetPlot {

  memory other_mem_req
  publishDir "${params.output_dir}/MissHetPlot", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
    path het 
    path imiss
  
  output:
    path output

  script:
    base = imiss.baseName
    output  = "${base}-imiss-vs-het".replace(".","_")+".pdf"
    template "missHetPlot.py"
}

// Find those who have bad heterozygosity
process getBadIndivsMissingHet {

  memory other_mem_req
  publishDir "${params.output_dir}/BadIndivsMissingHet", overwrite:true, mode:'copy', pattern: "*.txt"

  input:
    path het
    path imiss

  output:
    path outfname

  script:
    base = het.baseName
    outfname = "${base}-fail_het".replace(".","_")+".txt"
    template "select_miss_het_qcplink.py"
}

process removeQCIndivs {

  memory plink_mem_req

  input:
    path f_miss_het
    path rel_indivs

    path f_sex_check_f
    path sexcheck_report
    path hwe
    path imiss 
    path lmiss

    path poorgc
    path plates

    path bed
    path bim
    path fam
    path qc1
    path irem

  output:
    path "${out}.{bed,bim,fam}"

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

process calculateSnpSkewStatus {

  memory plink_mem_req
  cpus max_plink_cores

  input:
    path plinks

  output:
    path "${base}.missing"
    path mperm
    path "${base}.hwe"

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
   publishDir "${params.output_dir}/DMP", overwrite:true, mode:'copy', pattern: "*.pdf"

   input:
     path clean_missing
     path mperm
     path hwe
   
   output:
      path output

   script:
       input = clean_missing
       base  = clean_missing.baseName.replace(".","_").replace("-nd","")
       output= "${base}-diff-snpmiss_plot.pdf"
       template "diffMiss.py"
 }

process findSnpExtremeDifferentialMissingness {

  echo true
  memory other_mem_req

  input:
    path clean_missing_1
    path clean_missing
    path hwe

  output:
     path failed

  script:

    cut_diff_miss = params.cut_diff_miss
    missing = clean_missing
    base     = missing.baseName.replace("-.*","").replace(".","_")
    probcol = "EMP2"
    failed   = "sampleA-failed_diffmiss.snps"

    """echo $cut_diff_miss and $missing and $base and $probcol and $failed and 1000 """
    template "select_diffmiss_qcplink.py"
}

process removeSkewSnps {

  memory plink_mem_req
  publishDir "${params.output_dir}/skewSnps", overwrite:true, mode:'copy'  

  input:
    path plinks
    path failed
  
  output:
    path "${output}.bed"
    path "${output}.bim" 
    path "${output}.fam" 
    path "${output}.log"
  
  script:
  base = plinks[0].baseName
  output = params.output.replace(".","_")
  """
  plink $K --bfile $base $sexinfo --exclude $failed --make-bed --out $output
  """
}

process calculateMaf {

  memory plink_mem_req
  publishDir "${params.output_dir}/Maf", overwrite:true, mode:'copy', pattern: "*.frq"

  input:
    path bed
    path bim
    path fam
    path log

  output:
    path "${base}.frq"

  script:
    base = bed.baseName
    out  = base.replace(".","_")
    """
      plink --bfile $base $sexinfo  --freq --out $out
    """
}

process generateMafPlot {

  memory other_mem_req
  publishDir "${params.output_dir}/Maf", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
    path input

  output:
    path output

  script:
    base    = input.baseName
    output  = "${base}-maf_plot.pdf"
    template "mafplot.py"
}

process findHWEofSNPs {

  memory other_mem_req

  input:
    path missing
    path mperm
    path hwe

  output:
     path output

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
  publishDir "${params.output_dir}/HwePlot", overwrite:true, mode:'copy', pattern: "*.pdf"

  input:
    path unaff
  
  output:
    path output

  script:
    input  = unaff
    base   = unaff.baseName.replace(".","_")
    output = "${base}-hwe_plot.pdf"
    template "hweplot.py"
}

// Generate MD5 sums of output files
process outMD5 {
  echo true

  input:
     path bed
     path bim
     path fam
     path log

  output:
     path out

  script:
       out  = "${bed.baseName}.md5"
       template "md5.py"
}

process batchProc {

  memory plink_mem_req
  publishDir "${params.output_dir}/batchProc", pattern: "*{csv,pdf}", \
             overwrite:true, mode:'copy'

  input:
    path eigenval
    path eigenvec
    path bed
    path bim
    path fam

    path logfile
    path sexcheck_report
    path hwe
    path imiss
    path lmiss

    path "pheno.phe"
    path "batch.phe"
    path genome
    path pkl
    path rem_indivs

  
  output:
      path "${base}-batch.tex"
      path "*.csv"
      path "*pdf"

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
  publishDir "${params.output_dir}/batchProc", overwrite:true, mode:'copy'

  input:
    path bed
    path bim
    path fam
    path orig 
    path dupf
    path lmissf
    path imissf

    path cbed
    path cbim
    path cfam
    path ilog

    path missingvhetpdf

    path mafpdf

    path snpmisspdf

    path indmisspdf

    path fsex
    path sexcheck_report
    path hwe
  
    path misshetremf

    path diffmisspdf

    path diffmiss

    path eigenvalpdf
    path pcapdf

    path hwepdf

    path rel_indivs

    path inpmd5

    path outmd5

    path initmafpdf
    path initmaftex

    path inithwepdf
    path inithweqqpdf
    path inithwetex

    path bedd
    path bimm
    path famm
    path qc1
    path irem

    path batch_tex
    path bcsvs
    path bpdfs

    path poorgc10
    path poorgc

  output:
    path "${base}.pdf"

   script:
     base = "out"

     template "qcreport.py"
}