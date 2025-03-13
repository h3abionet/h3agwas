include {sampleSheet; getDuplicateMarkers; checkSampleSheet;clean_x;getX;identifyIndivDiscSexinfo;analyseX;;getInitMAF;showInitMAF;generateSnpMissingnessPlot;generateIndivMissingnessPlot; showHWEStats;removeQCPhase1;compPCA;drawPCA;pruneForIBDLD;findRelatedIndiv;calculateSampleHeterozygosity;generateMissHetPlot;getBadIndivsMissingHet;removeQCIndivs;calculateSnpSkewStatus;generateDifferentialMissingnessPlot;findSnpExtremeDifferentialMissingness;removeSkewSnps;computed_stat_female_x;computed_stat_male_x;report_export_x;cleanandmerge_x;report_export_x_tmp;clean_y;splitX;cleanPlink_x;cleanandmerge_y;build_reporty;calculateMaf;generateMafPlot;findHWEofSNPs;generateHwePlot;removeDuplicateSNPs;batchProc;produceReports} from './process.nf'
include {checkColumnHeader;fileheader_create_ch; strmem}  from '../modules/fct_groovy.nf'
include {is_nullfile; fileexist} from '../modules/fct_groovy.nf'
include {MD5_plk as MD5_in} from '../modules/utils.nf'
include {MD5_plk as MD5_out} from '../modules/utils.nf'
include {countChr as countX} from './process.nf'
include {countChr as countXY} from './process.nf'
include {countChr as countY} from './process.nf'


workflow check_params{
 take :
  data
  bfile 
 main : 
  filescript=file(workflow.scriptFile)                                            
  projectdir="${filescript.getParent()}" 
 /* check id*/
  if (workflow.repository)
    wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
  else
    wflowversion="A local copy of the workflow was used"
 if(params.pheno!='')pheno_col=params.pheno
 if(params.pheno_col!='')pheno_col=params.pheno_col
 if(params.pheno_qc!='')pheno_col=params.pheno_qc

 if(data){
  phenotype_ch = data
 }else {
   
   if(params.pheno!=""){
        phenotype=params.pheno
        println('warning :args phenotype will be deleted used --data or params.data')
   }else if (params.data!='')phenotype=params.data
   
   fileexist(phenotype)
   idfiles = [params.batch,phenotype]
   idfiles.each { checkColumnHeader(it,['FID','IID']) }
   phenotype_ch = fileheader_create_ch(phenotype,"pheno",pheno_col)
 }
 
 if(params.batch!='' && params.batch!='0')batch_ch     = fileheader_create_ch(params.batch,"batch",params.batch_col)
 else if(params.batch_col!=''){
   batch_ch=phenotype_ch
 }else{
  batch_ch = Channel.fromPath("$projectdir/assets/NO_FILE_4")
 }
 println(params.batch)
 diffpheno = ""
 if (params.case_control) {
  ccfile = params.case_control
  cc_ch = Channel.fromPath(ccfile,checkIfExists:true)
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
     if(params.case_control_col!=''){
       col    = params.case_control_col
       diffpheno = "--pheno cc.phe --pheno-name $col"
       cc_ch=phenotype_ch
     }else{
       col = ""
       cc_ch  = Channel.fromPath("$projectdir/assets/NO_FILE",checkIfExists:true)
    }
 }



 if (params.sexinfo_available) {
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

if(bfile){
  plink_ch=bfile
  bim_ch=bfile.flatMap{it -> it[1]}
}else {
 if(params.bfile!=''){
	inpat=params.bfile
 }else{
	inpat = "${params.input_dir}/${params.input_pat}"
 }
plink_ch=Channel.fromPath("${inpat}.bed",checkIfExists:true).combine(Channel.fromPath("${inpat}.bim",checkIfExists:true)).combine(Channel.fromPath("${inpat}.fam",checkIfExists:true))
bim_ch=Channel.fromPath("${inpat}.bim",checkIfExists:true)
}


if (params.idpat ==  "0")
    idpat   = "(.*)"
else
    idpat   = params.idpat

 if (is_nullfile(params.samplesheet)){
     samplesheet = "$projectdir/assets/NO_FILE_3"
     sample_sheet_ch = channel.fromPath(samplesheet,checkIfExists:true).combine(channel.of(0)).combine(channel.of(idpat))
 }else {
      samplesheet = params.samplesheet
   checkSampleSheet(samplesheet)
     sample_sheet_ch = channel.fromPath(samplesheet,checkIfExists:true).combine(channel.of(1)).combine(channel.of(idpat))
 }

 if(params.high_ld_regions_fname != "")   ldreg_ch=Channel.fromPath(params.high_ld_regions_fname,  checkIfExists:true) else ldreg_ch = Channel.fromPath("$projectDir/assets/NO_FILE",   checkIfExists:true)

  remove_on_bp    = channel.of(params.remove_on_bp)
  // todo check value
  pi_hat=channel.of(params.pi_hat)
 emit :
  wflowversion = wflowversion
  phenotype_ch = phenotype_ch
  pheno_col = pheno_col
  diffpheno = channel.of(diffpheno)
  batch_ch = batch_ch
  cc_ch = cc_ch
  col = col
  sexinfo = sexinfo
  extrasexinfo = extrasexinfo
  plink_ch = plink_ch
  bim_ch = bim_ch
  sample_sheet_ch= sample_sheet_ch
  remove_on_bp = remove_on_bp
  bal_sexinfo = bal_sexinfo
  pi_hat = pi_hat
ldreg_ch = ldreg_ch
}

workflow qc {
 take :                                                                         
  data                                                                          
  bfile                                                                         
  outputdir
 main :
 // check params and take param
 check_params(data,bfile)
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
 county=countY(removeDuplicateSNPs.out.plink, channel.of(params.chry_plink))
 // split x in 25 or 23  
 clean_x(removeDuplicateSNPs.out.plink, countX.out, countXY.out)
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
 identifyIndivDiscSexinfo(clean_x.out.plink)
 // update should be before or after removeDuplicateSNPs?
 generateSnpMissingnessPlot(removeDuplicateSNPs.out.lmiss)
 // update should be before or after removeDuplicateSNPs?
 generateIndivMissingnessPlot(removeDuplicateSNPs.out.imiss)
 //
 getInitMAF(clean_x.out.plink)
 showInitMAF(getInitMAF.out)
 showHWEStats(identifyIndivDiscSexinfo.out.hwe)
 removeQCPhase1(clean_x.out.plink, check_params.out.sexinfo)
 //
 compPCA(removeQCPhase1.out.plink)
 drawPCA(compPCA.out.eigen, check_params.out.cc_ch, check_params.out.col, check_params.out.diffpheno)
 pruneForIBDLD(removeQCPhase1.out.plink, check_params.out.ldreg_ch, check_params.out.sexinfo, check_params.out.pi_hat)
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
 clean_y(clean_x.out,removeQCIndivs.out.fam, county)
 cleanandmerge_y(dataf_withx_ch ,clean_y.out.plky_qc, county)
 build_reporty(clean_y.out.stat_qc, county)
 calculateMaf(removeSkewSnps.out, check_params.out.sexinfo) | generateMafPlot
 calculateSnpSkewStatus.out.hwe | findHWEofSNPs | generateHwePlot
 MD5_out(cleanandmerge_y.out.plk_log)
 batchProc(compPCA.out.eigen, identifyIndivDiscSexinfo.out.stat, check_params.out.phenotype_ch, check_params.out.batch_ch,pruneForIBDLD.out, x_analy_res_ch,findRelatedIndiv.out, check_params.out.extrasexinfo, check_params.out.pheno_col)
 produceReports(
      removeDuplicateSNPs.out.dups,
      cleanandmerge_y.out.plk_log,
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
 emit :
  plink=cleanandmerge_y.out.plk
}


include {dl_dataref_michigan;michigan_qc;clean_michigan} from './process_michigan.nf'
include {getfrequency} from '../modules/utils_plink.nf'
workflow qc_michigan {
  take :
   plink
  main :
   if(plink==null){
     if(params.bfile!=''){                                                          
        inpat=params.bfile                                                      
     }else{                                                                         
        inpat = "${params.input_dir}/${params.input_pat}"                       
      }                                                                              
   plink=Channel.fromPath("${inpat}.bed",checkIfExists:true).combine(Channel.fromPath("${inpat}.bim",checkIfExists:true)).combine(Channel.fromPath("${inpat}.fam",checkIfExists:true))
   }
   if (params.qc_michigan==1) {
      bin_checkmich=Channel.fromPath(params.bin_checkmich)
     if(params.dataref_michigan==""){
       dl_dataref_michigan()
       data_michigan =dl_dataref_michigan.out
     }else {
       data_michigan=  channel.fromPath(params.dataref_michigan, checkIfExists:true)
     }
    getfrequency(plink)
    michigan_qc(getfrequency.out, bin_checkmich, data_michigan)
    clean_michigan(michigan_qc.out.res, plink)
    plinkf = clean_michigan.out
   } else {
     plinkf=plink
  }
  emit :
    plink = plinkf
}
include {clean_phenofile;extract_ind_plink;computed_relatdness;plink_indep;check_rel} from './process_checkduplicate.nf'
include {computed_relatdness as computed_relatdness_dup} from './process_checkduplicate.nf'
include {computed_relatdness as computed_relatdness_all} from './process_checkduplicate.nf'
include {compute_missing as compute_missing_dup} from './process_checkduplicate.nf'
include {plink_updatename} from './process_checkduplicate.nf'

workflow qc_dup{
  take :
     data
     plink
     outputdir
  main :
   if(plink==null){
     if(params.bfile!=''){
        inpat=params.bfile
     }else{
        inpat = "${params.input_dir}/${params.input_pat}"
      }
     plink=Channel.fromPath("${inpat}.bed",checkIfExists:true).combine(Channel.fromPath("${inpat}.bim",checkIfExists:true)).combine(Channel.fromPath("${inpat}.fam",checkIfExists:true))
   } 
     
  if(data==null){
    data=channel.fromPath(params.data)
  }
  clean_phenofile(data, plink) 
  extract_ind_plink(plink,clean_phenofile.out.correspond)
  plink_indep(extract_ind_plink.out.combine(channel.of("${outputdir}/ind/")).combine(channel.of("${params.output}_indeprel")))
  computed_relatdness_dup(plink_indep.out.plk, channel.of(1).combine(clean_phenofile.out.dup), channel.of(0), channel.of(-1),channel.of("${params.output}_duplicate"))
  computed_relatdness_all(plink_indep.out.plk, channel.of(1).combine(clean_phenofile.out.dup), channel.of(0), channel.of(params.pi_hat_dup), channel.of("${params.output}_full"))
  compute_missing_dup(plink_indep.out.plk, channel.of(1).combine(clean_phenofile.out.dup),channel.of("${params.output}_duplicate"))
  check_rel(clean_phenofile.out.correspond, clean_phenofile.out.pheno_i, computed_relatdness_all.out, computed_relatdness_dup.out, compute_missing_dup.out.ind, channel.of(params.pi_hat_dup), channel.of("${params.output}_duplicate"), channel.of("$outputdir/checkrel"))
  plink_updatename(plink, check_rel.out.update_id, channel.of("${params.output}_qcdup"), channel.of("${outputdir}"))
  emit :
   data=check_rel.out.pheno
   plink=plink_updatename.out

}
