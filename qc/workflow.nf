include {sampleSheet; getDuplicateMarkers; checkSampleSheet;clean_x;getX;identifyIndivDiscSexinfo;analyseX;;getInitMAF;showInitMAF;generateSnpMissingnessPlot;generateIndivMissingnessPlot; showHWEStats;removeQCPhase1;compPCA;drawPCA;pruneForIBDLD;findRelatedIndiv;calculateSampleHeterozygosity;generateMissHetPlot;getBadIndivsMissingHet;removeQCIndivs;calculateSnpSkewStatus;generateDifferentialMissingnessPlot;findSnpExtremeDifferentialMissingness;removeSkewSnps;computed_stat_female_x;computed_stat_male_x;report_export_x;cleanandmerge_x;report_export_x_tmp;clean_y;splitX;cleanPlink_x;cleanandmerge_y;build_reporty;calculateMaf;generateMafPlot;findHWEofSNPs;generateHwePlot;removeDuplicateSNPs;batchProc;produceReports} from './process.nf'
include {checkColumnHeader;fileheader_create_ch; strmem}  from '../modules/fct_groovy.nf'
include {is_nullfile; fileexist} from '../modules/fct_groovy.nf'
include {MD5_plk as MD5_in} from '../modules/utils.nf'
include {MD5_plk as MD5_out} from '../modules/utils.nf'
include {countChr as countX} from './process.nf'
include {countChr as countXY} from './process.nf'
include {countChr as countY} from './process.nf'


workflow check_params{
 allowed_params= ["AMI","accessKey","batch","batch_col","bootStorageSize","case_control","case_control_col", "chipdescription", "cut_het_high","cut_get_low","cut_maf","cut_mind","cut_geno","cut_hwe","f_hi_female","f_lo_male","cut_diff_miss","cut_het_low", "help","input_dir","input_pat","instanceType","manifest", "maxInstances", "max_plink_cores","high_ld_regions_fname","low_memory","output", "output_align", "output_dir","phenotype","pheno_col","pi_hat", "plink_mem_req","region","reference","samplesheet", "scripts","secretKey","sexinfo_available", "sharedStorageMount","strandreport","work_dir","max_forks","big_time","super_pi_hat","samplesize","idpat","newpat","access-key","secret-key","instance-type","boot-storage-size","max-instances","shared-storage-mount","gemma_num_cores","remove_on_bp","queue","data","pheno","gc10", "build_genome", "autosome_plink", "cut_maf_xfemale", "cut_maf_xmale", "cut_miss_xfemale", "cut_miss_xmale", "cut_diffmiss_x", "chrxx_plink", "chrxy_plink", "chry_plink", "chrm_plink" ,"cut_maf_y", "cut_miss_y", "action", "bfile"]

 params.each { parm ->
  if (! allowed_params.contains(parm.key)) {
        println "Check $parm  ************** is it a valid parameter -- are you using one rather than two - signs or vice-versa";
   }
  }

  /* check id*/
 fileexist(params.batch)
 fileexist(params.phenotype)
 idfiles = [params.batch,params.phenotype]
 idfiles.each { checkColumnHeader(it,['FID','IID']) }

  if (workflow.repository)
    wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
  else
    wflowversion="A local copy of the workflow was used"

  phenotype_ch = fileheader_create_ch(params.phenotype,"pheno",params.pheno_col)
  batch_ch     = fileheader_create_ch(params.batch,"batch",params.batch_col)

 if (params.case_control) {
  ccfile = params.case_control
  cc_ch = Channel.fromPath(ccfile)
  col    = params.case_control_col
  diffpheno = "--pheno cc.phe --pheno-name $col"
  if (params.case_control.toString().contains("s3://") || params.case_control.toString().contains("az://")) {
       println "Case control file is in the cloud so we can't check it"
  } else
  if (! fileexist(params.case_control)) {
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
  cc_ch  = Channel.fromPath('00')
 }


 if (is_nullfile(params.sexinfo_available) ) {
   sexinfo = "--allow-no-sex"
   bal_sexinfo=false
   extrasexinfo = ""
   println "Sexinfo not available, command --allow-no-sex\n"
 } else {
   sexinfo = ""
   extrasexinfo = "--must-have-sex"
   println "Sexinfo available command"
   bal_sexinfo=true
 }

 if(params.bfile!=''){
	inpat=params.bfile
 }else{
	inpat = "${params.input_dir}/${params.input_pat}"

 }
plink_ch=Channel.fromPath("${inpat}.bed",checkIfExists:true).combine(Channel.fromPath("${inpat}.bim",checkIfExists:true)).combine(Channel.fromPath("${inpat}.fam",checkIfExists:true))
bim_ch=Channel.fromPath("${inpat}.bim",checkIfExists:true)


if (params.idpat ==  "0")
    idpat   = "(.*)"
else
    idpat   = params.idpat

 if (is_nullfile(params.samplesheet)){
     samplesheet = "0"
     sample_sheet_ch = channel.fromPath(samplesheet).combine(channel.of(0)).combine(channel.of(idpat))
 }else {
      samplesheet = params.samplesheet
drawPCA     checkSampleSheet(samplesheet)
     sample_sheet_ch = channel.fromPath(samplesheet).combine(channel.of(1)).combine(channel.of(idpat))
 }



  remove_on_bp    = channel.of(params.remove_on_bp)
  // todo check value
  pi_hat=channel.of(params.pi_hat)
 emit :
  wflowversion = wflowversion
  phenotype_ch = phenotype_ch   
  diffpheno = diffpheno
  batch_ch = batch_ch
  cc_ch = cc_ch
  col = col

  sexinfo = sexinfo
  extrasexinfo = extrasexinfo
  inpat = inpat
  plink_ch = plink_ch
  bim_ch = bim_ch
  sample_sheet_ch= sample_sheet_ch
  remove_on_bp = remove_on_bp
  bal_sexinfo = bal_sexinfo
  pi_hat = pi_hat
}

workflow qc {
 // check params and take param
 check_params()
 // compute mdmt
 MD5_in(check_params.out.plink_ch)
 sampleSheet(check_params.out.sample_sheet_ch)
 // defined using bim file duplicate marker
 getDuplicateMarkers(check_params.out.bim_ch, check_params.out.remove_on_bp)
 // remove duplicate marker
 removeDuplicateSNPs(check_params.out.plink_ch, getDuplicateMarkers.out, check_params.out.extrasexinfo, check_params.out.sexinfo)
 // count x number
 countxx=countX(removeDuplicateSNPs.out.plink, channel.of(params.chrxx_plink))
 // count xy number
 countxy=countXY(removeDuplicateSNPs.out.plink, channel.of(params.chrxy_plink))
 // count x number
 countY(removeDuplicateSNPs.out.plink, channel.of(params.chry_plink))
 clean_x(removeDuplicateSNPs.out.plink, countxy, countxx)
 cleanx=false
 if(check_params.out.bal_sexinfo){
  getX(clean_x.out) 
  analyseX(clean_x.out)
  cleanx=true
   x_analy_res_ch=analyseX.out
 }else {
   //need to update
   x_analy_res_ch = Channel.fromPath("0")
 }
 identifyIndivDiscSexinfo(removeDuplicateSNPs.out.plink)
 // update should be before or after removeDuplicateSNPs? 
 generateSnpMissingnessPlot(removeDuplicateSNPs.out.lmiss)
 // update should be before or after removeDuplicateSNPs? 
 generateIndivMissingnessPlot(removeDuplicateSNPs.out.imiss)
 //
 getInitMAF(removeDuplicateSNPs.out.plink)
 showInitMAF(getInitMAF.out)
 showHWEStats(identifyIndivDiscSexinfo.out.hwe)
 removeQCPhase1(removeDuplicateSNPs.out.plink, check_params.out.sexinfo)
 // 
 compPCA(removeQCPhase1.out.plink)
 drawPCA(compPCA.out.eigen, check_params.out.cc_ch, check_params.out.col, check_params.out.diffpheno)
 if(params.high_ld_regions_fname != "")   ldreg_ch=Channel.fromPath(params.high_ld_regions_fname,  checkIfExists:true) else ldreg_ch = Channel.fromPath("$projectDir/assets/NO_FILE",   checkIfExists:true)
 pruneForIBDLD(removeQCPhase1.out.plink, ldreg_ch, check_params.out.sexinfo, check_params.out.pi_hat)
 // why function used  previous  remove duplicate miss?
 findRelatedIndiv(removeDuplicateSNPs.out.imiss,pruneForIBDLD.out)
 calculateSampleHeterozygosity(removeQCPhase1.out.plink, check_params.out.sexinfo)
 generateMissHetPlot(calculateSampleHeterozygosity.out.hetmiss)
 getBadIndivsMissingHet(calculateSampleHeterozygosity.out.hetmiss)
 removeQCIndivs(getBadIndivsMissingHet.out, findRelatedIndiv.out, identifyIndivDiscSexinfo.out.log,sampleSheet.out.poorsgc10,removeQCPhase1.out.plink, check_params.out.sexinfo)
 calculateSnpSkewStatus(removeQCIndivs.out.plink,check_params.out.cc_ch, check_params.out.sexinfo, check_params.out.diffpheno)
 generateDifferentialMissingnessPlot(calculateSnpSkewStatus.out.missing)
 findSnpExtremeDifferentialMissingness(calculateSnpSkewStatus.out.perm)
 removeSkewSnps(removeQCIndivs.out.plink, findSnpExtremeDifferentialMissingness.out.failed, check_params.out.sexinfo)
 if(cleanx){
  splitX(clean_x.out,removeQCIndivs.out.fam)
  cleanPlink_x(splitX.out)
  computed_stat_female_x(cleanPlink_x.out.female_x)
  computed_stat_male_x(cleanPlink_x.out.male_x)
  report_export_x(computed_stat_male_x.out, computed_stat_female_x.out)
  cleanandmerge_x(removeSkewSnps.out, splitX.out,report_export_x.out.listsnps)
  report_x = report_export_x.out.report
  dataf_withx_ch = cleanandmerge_x.out
 }else{
  report_export_x_tmp()
  report_x = report_export_x_tmp.out 
  dataf_withx_ch = removeSkewSnps.out
 }
  
 clean_y(splitX.out,removeQCIndivs.out.fam, countY.out)  
 cleanandmerge_y(dataf_withx_ch ,clean_y.out.plky_qc, countY.out)
 build_reporty(clean_y.out.stat_qc, countY.out) 
 calculateMaf(removeSkewSnps.out, check_params.out.sexinfo) | generateMafPlot
 calculateSnpSkewStatus.out.hwe | findHWEofSNPs | generateHwePlot
 MD5_out(cleanandmerge_y.out)
 batchProc(compPCA.out.eigen, identifyIndivDiscSexinfo.out.stat, check_params.out.phenotype_ch, check_params.out.batch_ch,pruneForIBDLD.out, x_analy_res_ch,findRelatedIndiv.out, check_params.out.extrasexinfo)
 produceReports(
      removeDuplicateSNPs.out.dups, 
      cleanandmerge_y.out, 
      generateMissHetPlot.out, 
      generateMafPlot.out, 
      generateSnpMissingnessPlot.out, 
      generateIndivMissingnessPlot.out, 
      identifyIndivDiscSexinfo.out.log,
      getBadIndivsMissingHet.out,
      generateDifferentialMissingnessPlot.out,
      findSnpExtremeDifferentialMissingness.out.failed,
      drawPCA.out,
      generateHwePlot.out,
      findRelatedIndiv.out, 
      MD5_in.out, 
      MD5_out.out, 
      showInitMAF.out,
      showHWEStats.out,
      removeQCPhase1.out.log,
      batchProc.out.report_batch_report_ch,
      sampleSheet.out.plates, 
      batchProc.out.report_batch_aux_ch,
      report_x, 
      build_reporty.out.report)
}
