params.diff_pheno=""
params.case_control_col=""
params.pheno="" // deprecated

pp = [:]



include {
   analyseX;
   batchProc;
   build_reporty;
   calculateMaf;
   calculateSampleHeterozygosity;
   calculateSnpSkewStatus;
   checkSampleSheet;
   cleanPlink_x;
   clean_y;
   cleanandmerge_x;
   cleanandmerge_y;
   compPCA;
   computed_stat_female_x;
   computed_stat_male_x;
   countChr as countX;
   countChr as countXY;
   countChr as countY
   drawPCA;
   findHWEofSNPs;
   findRelatedIndiv;
   findSnpExtremeDifferentialMissingness;
   generateDifferentialMissingnessPlot;
   generateHwePlot;
   generateIndivMissingnessPlot;
   generateMafPlot;
   generateMissHetPlot;
   generateSnpMissingnessPlot;
   getBadIndivsMissingHet;
   getDuplicateMarkers;
   getInitMAF;
   getX;
   identifyIndivDiscSexinfo;
   produceReports
   pruneForIBDLD;
   relabel_par;
   removeDuplicateSNPs;
   removeQCIndivs;
   removeQCPhase1;
   removeSkewSnps;
   report_export_x;
   report_export_x_tmp;
   sampleSheet;
   showInitMAF;
   showHWEStats;
   splitX;
} from './qc-processes.nf'


include {
    check_rel;
    clean_phenofile;
    compute_missing as compute_missing_dup;
    computed_relatdness as computed_relatdness_dup;
    computed_relatdness as computed_relatdness_all;
    extract_ind_plink;computed_relatdness;
    plink_indep;
    plink_updatename;
} from './check-duplicate-processes.nf'


include {checkColumnHeader;
	 fileheader_create_ch;
	 fileexist;
	 is_nullfile;
	 strmem}  from '../modules/fct_groovy.nf'


include {MD5_plk as MD5_in; MD5_plk as MD5_out} from '../modules/utils.nf'


// Code to check parameters


def try_cc_file(ccfile,cc_col) {
    cc_ch = Channel.fromPath(ccfile,checkIfExists:true)
    if (params.case_control.toString().contains("s3://") ||
        params.case_control.toString().contains("az://")) {
          println "Case control file is in the cloud so we can't check it"
	  return true
    }
    def line
    new File(ccfile).withReader { line = it.readLine() }
    fields = line.split()
    if (! fields.contains(cc_col)) {
	error("\n\nThe file <$ccfile}> given for <params.case_control>"+
	      "does not have a column"+"<${cc_col}>\n")
	return false
    } else {
	return true
    }
}
              

def check_case_control_params(pp) {
    diff_pheno = ""
    col        = params.case_control_col
    if (col == '') return  [diff_pheno, Channel.fromPath("$projectdir/assets/NO_FILE",checkIfExists:true)]
    diff_pheno = "--pheno cc.phe --pheno-name $col"
    cc = params.data
    if (col && try_cc_file(cc, col)) {
	pp['diff_pheno']=diff_pheno
	pp['cc_ch']=Channel.fromPath(cc)
       return;
    }
    error "Case control col <$col> is given but not found in <$cc> or $cc does not exist"
    System.exit(10)
}

def check_phenotype_file(data,pheno_col) {
    if (data){  // if so, then this is called as sub-sub-worfklow
        return data
    }
    if(params.pheno!="") {
	phenotype=params.pheno
	error('warning :args phenotype has  been deleted used --data or params.data')
	System.exit(12)
    } 
    if (params.data='') {
	error "No phenotype file specified"
	System.exit(13)
    }
    phenotype=params.data
    fileexist(phenotype)
    idfiles = [params.batch,phenotype]
    idfiles.each { checkColumnHeader(it,['FID','IID']) }
    return fileheader_create_ch(phenotype,"pheno",pheno_col)
}

def check_sex_params (pp) {
   if (!params.sex_info_available) {
	pp['sex_info'] = "--allow-no-sex"
        pp['bal_sex_info']=false
        pp['extra_sex_info'] = ""
        println "Sex_Info not available, command --allow-no-sex\n"
   } else {
        pp['sex_info'] = ""
        pp['extra_sex_info'] = "--must-have-sex"
        println "Sex_Info available command"
        pp['bal_sex_info']=true
   }
}

def check_plink_params(bfile,pp) {
   if(bfile){
	pp['plink_ch'] = bfile
        pp['bim_ch']   = bfile.flatMap{it -> it[1]}
   } else {
        if(params.bfile!=''){
	   inpat=params.bfile
        } else {
	   inpat = "${params.input_dir}/${params.input_pat}"
        }
	println inpat
        pp['plink_ch'] =  Channel.fromFilePairs("${inpat}.{bed,bim,fam}",size:3) { file -> file.baseName }\
                .ifEmpty { error "No matching ${inpat} bed,bim,fam file" }\
	        .map { it -> it[1] }
        pp['bim_ch']   =  Channel.fromPath("${inpat}.bim",checkIfExists:true)
    }
    
}



def check_params(data, bfile) {
    pp['filescript']=file(workflow.scriptFile)                                            
    pp['projectdir']="${pp['filescript'].getParent()}"
    pp['wflowversion']=  (workflow.repository) ? "${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
                                                     : "A local copy of the workflow was used"
    pp['pheno']     = params.pheno ?: null
    pp['pheno_col'] = params.pheno_col ?: null
    pp['pheno_qc']  = params.pheno_qc ?: null

    pp['phenotype_ch']=check_phenotype_file(data, pp['pheno_col'])
 
   if(params.data!='' && params.data!='0' && params.batch_col)
     batch_ch = fileheader_create_ch(params.data,"batch",params.batch_col)
    else
     batch_ch = Channel.fromPath("$projectdir/assets/NO_FILE_4")
   check_case_control_params(pp) 
   check_sex_params(pp)
   check_plink_params(bfile,pp)
   pp['id_pat'] = (params.idpat ==  "0") ? "(.*)" : params.idpat
   if  (is_nullfile(params.samplesheet)){
	samplesheet = "${pp['projectdir']}/assets/NO_FILE_3"
	pp['sample_sheet_ch'] = channel.fromPath(samplesheet,checkIfExists:true)
   } else  {
	samplesheet = params.samplesheet
	checkSampleSheet(samplesheet)
	pp['sample_sheet_ch'] = channel.fromPath(samplesheet,checkIfExists:true)
   }
   if(params.high_ld_regions_fname != "")
	pp['ldreg_ch']=Channel.fromPath(params.high_ld_regions_fname,  checkIfExists:true)
   else
	pp['ldreg_ch'] = Channel.fromPath("$projectDir/assets/NO_FILE",   checkIfExists:true)
   pp['diff_pheno'] = diff_pheno
}

workflow qc {
   take :                                                                         
    data                                                                          
    bfile                                                                         
    outputdir
   main :
    check_params(data,bfile) // sets pp
    // compute mdmt
    MD5_in(pp['plink_ch'])

    sampleSheet(pp['id_pat'],pp['sample_sheet_ch'])
   // defined using bim file duplicate marker
    getDuplicateMarkers(pp['bim_ch'])
    // remove duplicate marker
    removeDuplicateSNPs(pp['plink_ch'], getDuplicateMarkers.out,extra_sex_info:pp['extra_sex_info'],sex_info:pp['sex_info'])


   // count x number
    countxx=countX(params.chrxx_plink,removeDuplicateSNPs.out.plink)
    // count xy number
    countxy=countXY(params.chrxy_plink,removeDuplicateSNPs.out.plink)
    // count x number
    county=countY(params.chry_plink,removeDuplicateSNPs.out.plink)
    // split x in 25 or 23  

    
    relabelled = relabel_par(removeDuplicateSNPs.out.plink, countX.out, countXY.out)

    relabelled | identifyIndivDiscSexinfo
    relabelled | getInitMAF | showInitMAF
    relabelled | combine([pp['sex_info']]) |  removeQCPhase1 
    removeQCPhase1.out.plink | compPCA
    drawPCA(compPCA.out.eigen,  pp['cc_ch'], params.case_control_col)    

    x_analyse_res_ch = (params.sex_info_available) ? relabelled | (analyseX & getX) : Channel.fromPath("0")

    x_analyse_res_ch = analyseX.out

    generateSnpMissingnessPlot(removeDuplicateSNPs.out.lmiss)
    generateIndivMissingnessPlot(removeDuplicateSNPs.out.imiss)
    showHWEStats(identifyIndivDiscSexinfo.out.hwe)


    pruneForIBDLD(removeQCPhase1.out.plink, pp['ldreg_ch'], pp['sex_info'], params.pi_hat)

    related_indivs = findRelatedIndiv(removeDuplicateSNPs.out.imiss,pruneForIBDLD.out)

    calculateSampleHeterozygosity(removeQCPhase1.out.plink, pp["sex_info"]).hetmiss |  (generateMissHetPlot & getBadIndivsMissingHet)

    removeQCIndivs(getBadIndivsMissingHet.out, related_indivs, identifyIndivDiscSexinfo.out.log,
		   sampleSheet.out.poorsgc10, removeQCPhase1.out.plink, pp["sex_info"])

    calculateSnpSkewStatus(removeQCIndivs.out.plink,pp['cc_ch'], pp['sex_info'], pp['diff_pheno'])
    generateDifferentialMissingnessPlot(calculateSnpSkewStatus.out.missing)
    findSnpExtremeDifferentialMissingness(calculateSnpSkewStatus.out.perm)
    removeSkewSnps(removeQCIndivs.out.plink, findSnpExtremeDifferentialMissingness.out.failed, pp['sex_info'])
    if(params.sex_info_available){
	splitX(relabel_par.out,removeQCIndivs.out.fam)
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
   clean_y(relabel_par.out,removeQCIndivs.out.fam, county)
   cleanandmerge_y(dataf_withx_ch ,clean_y.out.plky_qc, county)
   build_reporty(clean_y.out.stat_qc, county)
   calculateMaf(removeSkewSnps.out, pp['sex_info'])| generateMafPlot
   calculateSnpSkewStatus.out.hwe | findHWEofSNPs | generateHwePlot
   MD5_out(cleanandmerge_y.out.plk_log)

   batchProc(compPCA.out.eigen, identifyIndivDiscSexinfo.out.stat, pp['phenotype_ch'], batch_ch,\
	     pruneForIBDLD.out, x_analyse_res_ch, findRelatedIndiv.out, pp['extra_sex_info'],pp['pheno_col'])

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
