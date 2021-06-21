#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {

    inMD5;
    sampleSheet;
    noSampleSheet; 
    getDuplicateMarkers;
    removeDuplicateSNPs;
    getX;
    analyseX;
    identifyIndivDiscSexinfo;
    identifyIndivDiscSexinfo_2;
    generateSnpMissingnessPlot;
    generateIndivMissingnessPlot;
    getInitMAF;
    showInitMAF;
    showHWEStats;
    removeQCPhase1;
    compPCA;
    drawPCA;
    pruneForIBD;
    findRelatedIndiv;
    calculateSampleHeterozygosity;
    generateMissHetPlot;
    getBadIndivsMissingHet;
    removeQCIndivs;
    calculateSnpSkewStatus;
    generateDifferentialMissingnessPlot;
    findSnpExtremeDifferentialMissingness;
    removeSkewSnps;
    calculateMaf;
    generateMafPlot;
    findHWEofSNPs;
    generateHwePlot;
    outMD5;
    batchProc;
    produceReports;

} from './modules/MULTIQC.nf';

// if (!workflow.resume) {
//     def dir = new File(params.output_dir)
//     if (dir.exists() && dir.directory && (!(dir.list() as List).empty)) {
//        println "\n\n============================================"
//        println "Unless you are doing a -resume, the output directory should be empty"
//        println "We do not want to overwrite something valuable in "+params.output_dir
//        println "Either clean your output directory or check if you meant to do a -resume"
//        System.exit(-1)
//     }
// }

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

nextflowversion =nextflow.version

if (workflow.repository)
  wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
else
  wflowversion="A local copy of the workflow was used"

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


if ( nullfile.contains(params.sexinfo_available) ) {
  sexinfo = "--allow-no-sex"
  extrasexinfo = ""
} else {
  sexinfo = ""
  extrasexinfo = "--must-have-sex"
}

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






workflow {

    inMD5(checked_input.inpmd5ch)

    println samplesheet
    if (samplesheet != "0")  {

        sample_sheet_ch = file(samplesheet)

        sampleSheet(sample_sheet_ch)

    }else{

        noSampleSheet()
    }

    getDuplicateMarkers(checked_input.bim_ch)

    removeDuplicateSNPs(checked_input.raw_ch,getDuplicateMarkers.out)

    missingness = [0.01,0.03,0.05]

    if (extrasexinfo == "--must-have-sex") {
        getX(removeDuplicateSNPs.out)
        analyseX(getX.out)
    }else{
        analyseX.out = Channel.fromPath("0")
    }

    identifyIndivDiscSexinfo(removeDuplicateSNPs.out)

    identifyIndivDiscSexinfo_2(removeDuplicateSNPs.out)

    generateSnpMissingnessPlot(removeDuplicateSNPs.out)

    generateIndivMissingnessPlot(removeDuplicateSNPs.out)

    getInitMAF(removeDuplicateSNPs.out)

    showInitMAF(getInitMAF.out)

    showHWEStats(identifyIndivDiscSexinfo.out)

    removeQCPhase1(removeDuplicateSNPs.out)

    compPCA(removeQCPhase1.out)

    drawPCA(compPCA.out, cc_ch)

    pruneForIBD(removeQCPhase1.out)

    findRelatedIndiv(removeDuplicateSNPs.out ,pruneForIBD.out)

    calculateSampleHeterozygosity(removeQCPhase1.out)

    generateMissHetPlot(calculateSampleHeterozygosity.out)

    getBadIndivsMissingHet(calculateSampleHeterozygosity.out)

    if (samplesheet != "0")  {

        removeQCIndivs(getBadIndivsMissingHet.out,findRelatedIndiv.out
    ,identifyIndivDiscSexinfo.out,sampleSheet.out, removeQCPhase1.out )

    }else{

        removeQCIndivs(getBadIndivsMissingHet.out,findRelatedIndiv.out
    ,identifyIndivDiscSexinfo.out,noSampleSheet.out, removeQCPhase1.out )

    }

    calculateSnpSkewStatus(removeQCIndivs.out.combine(cc_ch))

    generateDifferentialMissingnessPlot(calculateSnpSkewStatus.out)

    findSnpExtremeDifferentialMissingness(calculateSnpSkewStatus.out)

    removeSkewSnps(removeQCIndivs.out,findSnpExtremeDifferentialMissingness.out)

    calculateMaf(removeSkewSnps.out)

    generateMafPlot(calculateMaf.out)

    findHWEofSNPs(calculateSnpSkewStatus.out)

    generateHwePlot(findHWEofSNPs.out)

    outMD5(removeSkewSnps.out)

    batchProc(compPCA.out,identifyIndivDiscSexinfo.out
    ,phenotype_ch,batch_ch,pruneForIBD.out,analyseX.out,findRelatedIndiv.out)

    if (samplesheet != "0")  {

      produceReports(removeDuplicateSNPs.out,removeSkewSnps.out,generateMissHetPlot.out,
      generateMafPlot.out,generateSnpMissingnessPlot.out,generateIndivMissingnessPlot.out,
      identifyIndivDiscSexinfo_2.out,getBadIndivsMissingHet.out,
      generateDifferentialMissingnessPlot.out,findSnpExtremeDifferentialMissingness.out,
      drawPCA.out,generateHwePlot.out,findRelatedIndiv.out,inMD5.out,outMD5.out,
      showInitMAF.out,showHWEStats.out,removeQCPhase1.out,batchProc.out,sampleSheet.out)

    }else{

      produceReports(removeDuplicateSNPs.out,removeSkewSnps.out,generateMissHetPlot.out,
      generateMafPlot.out,generateSnpMissingnessPlot.out,generateIndivMissingnessPlot.out,
      identifyIndivDiscSexinfo_2.out,getBadIndivsMissingHet.out,
      generateDifferentialMissingnessPlot.out,findSnpExtremeDifferentialMissingness.out,
      drawPCA.out,generateHwePlot.out,findRelatedIndiv.out,inMD5.out,outMD5.out,
      showInitMAF.out,showHWEStats.out,removeQCPhase1.out,batchProc.out,noSampleSheet.out)

    }

    final_ch = produceReports.out;

    final_ch.subscribe { b=it.baseName; 

    new File("emptyZ0batch.txt").delete();
    new File("emptyZ0pheno.txt").delete();
    new File("xxemptyZ0pheno.txt").delete();

    println "The output report is called ${params.output_dir}/batchProc/${b}.pdf"}
}