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
  main :
   if(vcf==null){
      if(params.vcf_folder_zip!=''){
         unzip_folder(channel.fromPath(params.vcf_folder_zip,checkIfExists:true))
         vcf=unzip_folder.out.vcf
     }else {
        if(params.vcf!='')vcf=channel.fromPath(params.vcf)
        if(params.vcf_list!='')vcf=Channel.fromPath(file(params.vcf_list).readLines(), checkIfExists:true)
     }
  }
  dl_fasta_wf(params.fasta, build, params.ftp_fasta, '')
 emit :
   vcf=vcf
   fasta_index = dl_fasta_wf.out.fasta_index
   fasta = dl_fasta_wf.out.fasta

}
workflow convertvcfinplk {
  take :
     vcf
     build
     outputdir
     outputpat
  main :
      convert_inplink(vcf.combine(channel.of("${outputpat}")).combine(channel.of("${outputdir}/plink")))
      GetRsDup(convert_inplink.out.bim.collect())
      TransformRsDup(GetRsDup.out.combine(convert_inplink.out.bed))
      plink=TransformRsDup.out.plk
     if(params.genetic_maps!=""){
      AddedCM(plink)
      plink=AddedCM.out.plk
     }
     convert_inplink.out.bim.view()
     /*if(convert_inplink.out.bim.count().toInteger()>1) {
        MergePlink(plink.collect(), channel.of("${outputdir}/plink"), channel.of("${outputpat}_clean"))
        plink = MergePlink.out
     }*/
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
  main :
    getparams(vcf, build)
    computedstat(getparams.out.vcf)
    dostat(computedstat.out.collect(), channel.of(outputdir+'/stats/'), channel.of(outputpat))
    latex_compilation(dostat.out.tex, dostat.out.support_file, channel.of(outputdir+'/report/'))
    clean_vcf(getparams.out.vcf.combine(getparams.out.fasta).combine(channel.of(params.cut_maf)).combine(channel.of(params.cut_hwe)).combine(channel.of(params.vcf_minscoreimp)).combine(channel.of(params.vcf_cut_miss)).combine(channel.of("$outputdir/cleanvcf")))
    if(params.convertvcfinplink==1) {
     convertvcfinplk(clean_vcf.out,build,"$outputdir/plink/",outputpat)
    }
 }


