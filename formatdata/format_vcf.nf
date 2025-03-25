include {unzip_folder;computedstat;dostat;MergePlink} from './process_vcf_format.nf'
include {bcftools_index_vcf} from '../modules/process_vcf.nf'
include {dl_fasta_wf} from '../modules/dl_wf.nf'
include {latex_compilation} from '../modules/utils.nf'
include {AddedCM;GetRsDup;TransformRsDup;check_names_plkconvert} from './process_vcf_format.nf'
include {convert_inplink} from './process_vcf_format.nf'
include {clean_vcf} from './process_vcf_format.nf'
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
     println('using vcf from vcf '+params.vcf+' or list_vcf '+params.vcf_list)
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
 if(vcf_imputeformat==null)vcf_imputeformat=params.vcf_imputeformat
 else vcf_imputeformat=vcf_imputeformat.val
 if(vcf_imputeformat.toLowerCase()=="pbwt"){
  println("[FORMATVCFINPLINK] using pbwt")
  vcf_patscoreimp= "INFO"
  vcf_patstatfreq="%AN %AC"
 }
 else if(vcf_imputeformat.toLowerCase()=="minmac4"){
  println("[FORMATVCFINPLINK] using minmac4")
  vcf_patscoreimp= "R2"
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
 bcftools_index_vcf(vcf)
 emit :
   vcf= vcf
   vcf_index=bcftools_index_vcf.out
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
    vcf_patscoreimp
  main :
      convert_inplink(vcf ,output_pat:"", output_dir:"${outputdir}/plink", patscoreimp:vcf_patscoreimp)
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
    println('[FORMATVCFINPLINK] convert vcf in ')
    if(outputdir==null){
       outputdir=params.output_dir
    }
    if(outputpat==null){
       outputpat = params.output
    }
    getparams(vcf, build,vcf_imputeformat)
    if(getparams.out.do_stat.val==1){
     println('[FORMATVCFINPLINK] perform stqatistics')
     computedstat(getparams.out.vcf, getparams.out.vcf_patstatfreq, getparams.out.vcf_patscoreimp)
     dostat(computedstat.out.collect(), output_dir:"$outputdir/stats/", output_pat:outputpat)
     latex_compilation(dostat.out.tex, dostat.out.support_file, channel.of(outputdir+'/report/'))
     }
    clean_vcf(getparams.out.vcf_index.combine(getparams.out.fasta_index),maf:params.cut_maf_postimp, hwe:params.cut_hwe, r2lim:params.impute_info_cutoff,miss:params.vcf_cut_miss,output_dir:"$outputdir/cleanvcf",patscoreimp:getparams.out.vcf_patscoreimp.val)
    if(params.convertvcfinplink==1) {
    convertvcfinplk(clean_vcf.out,build,"$outputdir/plink/",outputpat, getparams.out.vcf_patscoreimp.val)
    }
    emit :
     plink = convertvcfinplk.out.plink 
     vcf = clean_vcf.out 
 }


