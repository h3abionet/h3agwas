include {unzip_folder;computedstat;dostat;convert_inplink;MergePlink} from './process_vcf_format.nf'
include {dl_fasta_wf} from '../modules/dl_wf.nf'
include {latex_compilation} from '../modules/utils.nf'
//formatdata/process_vcf_format.nf
include {clean_vcf;AddedCM;GetRsDup;TransformRsDup;check_names_plkconvert} from './process_vcf_format.nf'
include {updateplk_rsname} from './process.nf'

workflow getparams{
  take :
   vcf
   build
   vcf_imputeformat
  main :
   if(vcf){
     println('using vcf from previous process')
     vcf=vcf
   }else {
     println('using vcf from vcf '+params.vcf+' or list_vcf '+params.list_vcf)
      if(params.vcf!='')vcf=channel.fromPath(params.vcf, checkIfExists:true)
      if(params.vcf_list!='')vcf=Channel.fromPath(file(params.vcf_list).readLines(), checkIfExists:true)
      if(params.vcf_folder_zip!=''){
         unzip_folder(channel.fromPath(params.vcf_folder_zip,checkIfExists:true))
         vcf=unzip_folder.out.vcf
     }
  }
   fasta=params.fasta
   if(params.fasta_convert!='')fasta=params.fasta_convert
  dl_fasta_wf(fasta, build, params.ftp_fasta, '')
 vcf_patscoreimp =""
 vcf_patstatfreq =""
 if(vcf_imputeformat.val.toLowerCase()=="pbwt"){
  println("[FORMATVCFINPLINK] using pbwt")
  vcf_patscoreimp= "%INFO/R2"
  vcf_patstatfreq="%AN %AC"
 }
 else if(vcf_imputeformat.val.toLowerCase()=="minmac4"){
  println("[FORMATVCFINPLINK] using minmac4")
  vcf_patscoreimp= "%INFO/R2"
  vcf_patstatfreq="%MAF"
 } 

 if(params.vcf_patscoreimp!=''){
   vcf_patscoreimp = params.vcf_patscoreimp
 }
 if(params.vcf_patstatfreq!=''){
   vcf_patstatfreq = params.vcf_patstatfreq
 }
 do_stat=0
 if(vcf_patscoreimp!='' && vcf_patstatfreq!='' && params.convertvcf_stat==1){
    println('performing statistics using vcf') 
   do_stat=1
 }
 emit :
   vcf=vcf
   fasta_index = dl_fasta_wf.out.fasta_index
   fasta = dl_fasta_wf.out.fasta
   vcf_patscoreimp = vcf_patscoreimp
   vcf_patstatfreq = vcf_patstatfreq
   do_stat = do_stat

}
workflow convertvcfinplk {
  take :
     vcf
     build
     outputdir
     outputpat
  main :
      convert_inplink(vcf.combine(channel.of("$outputpat")).combine(channel.of("${outputdir}/plink")).combine(vcf.count()))
      GetRsDup(convert_inplink.out.bim.collect())
      TransformRsDup(GetRsDup.out.combine(convert_inplink.out.bed))
      plink=TransformRsDup.out.plk
     if(params.genetic_maps!=""){
      AddedCM(plink)
      plink=AddedCM.out.plk
     }
    /**/
    MergePlink(plink.collect(), channel.of("${outputdir}/plink"), channel.of("${outputpat}_clean"))
    plink = MergePlink.out
     if(params.data!=""){
        check_names_plkconvert(plink, channel.fromPath(params.data, checkIfExists:true))
       plink = check_names_plkconvert.out
    }
  emit :
     plink = plink
}


workflow convertvcfin{
  take :
     vcf
     build
     outputdir
     outputpat
     vcf_imputeformat
  main :
   vcf_imputeformat.view()
    getparams(vcf, build,vcf_imputeformat)
    if(getparams.out.do_stat.val==1){
     computedstat(getparams.out.vcf, getparams.out.vcf_patstatfreq, getparams.out.vcf_patscoreimp)
     dostat(computedstat.out.collect(), channel.of(outputdir+'/stats/'), channel.of(outputpat))
     latex_compilation(dostat.out.tex, dostat.out.support_file, channel.of(outputdir+'/report/'))
     }
    clean_vcf(getparams.out.vcf.combine(getparams.out.fasta).combine(channel.of(params.cut_maf)).combine(channel.of(params.cut_hwe)).combine(channel.of(params.impute_info_cutoff)).combine(channel.of(params.vcf_cut_miss)).combine(channel.of("$outputdir/cleanvcf")))
    if(params.convertvcfinplink==1) {
    convertvcfinplk(clean_vcf.out,build,"$outputdir/plink/",outputpat)
    }
    emit :
     plink = convertvcfinplk.out.plink 
     vcf = clean_vcf.out 
 }


