include {bcftools_index_vcf} from '../modules/process_vcf.nf'
include {list_chro} from '../modules/utils_plink.nf'
include {updateplk_rsname;updateplk_rsname_norsinfo;deletedmultianddel;refallele;convertInVcf;checkfixref;checkVCF;mergevcf;get_sex_forsanger} from './process.nf'
include {fileexist_b;fileexist_param;fileexist} from '../modules/fct_groovy.nf'
include {dl_fasta_wf} from '../modules/dl_wf.nf'

/*workflow getparams{
  take :
   plink
   fasta
  main :
   if(!plink){
     if(params.bfile!=''){
        inpat=params.bfile
     }else{
        inpat = "${params.input_dir}/${params.input_pat}"
      }
   plink=Channel.fromPath(fileexist("${inpat}.bed"),checkIfExists:true).combine(Channel.fromPath(fileexist("${inpat}.bim"),checkIfExists:true)).combine(Channel.fromPath(fileexist("${inpat}.fam"),checkIfExists:true))
   }
 rs_infogz=null
 if(fileexist_b(params.file_ref_gzip)){
   if(fileexist_param(params.file_ref_gzip, "file_ref_gzip"))file_ref_gzip=channel.fromPath(params.file_ref_gzip)
     if(fileexist_b("${params.file_ref_gzip}.csi")){
     rs_infogz=file_ref_gzip.combine(channel.fromPath(params.file_ref_gzip+'.csi',checkIfExists:true))
   }else {
     bcftools_index_vcf(file_ref_gzip)
     rs_infogz=bcftools_index_vcf.out
   }
 }
 dl_fasta_wf(params.fasta, params.build_genome, params.ftp_fasta, '')

 emit :
  rs_infogz= rs_infogz
  fasta = dl_fasta_wf.out.fasta_index
  plink = plink
}
*/
pp=[:]
def check_params(bfile, fasta) {
 pp['bfile'] = bfile
 if (bfile==null){
     if(params.bfile!=''){
        inpat=params.bfile
     }else{
        inpat = "${params.input_dir}/${params.input_pat}"
      }
   pp['bfile'] = Channel.fromPath(fileexist("${inpat}.bed"),checkIfExists:true).combine(Channel.fromPath(fileexist("${inpat}.bim"),checkIfExists:true)).combine(Channel.fromPath(fileexist("${inpat}.fam"),checkIfExists:true))
 }
 pp['fasta'] = fasta
 if(fasta == null){
   dl_fasta_wf(params.fasta, params.build_genome, params.ftp_fasta, '')
  pp['fasta'] =  dl_fasta_wf.out.fasta_index
 }
 pp['rs_info'] = null
 if(fileexist_b(params.file_ref_gzip)){
   if(fileexist_param(params.file_ref_gzip, "file_ref_gzip"))file_ref_gzip=channel.fromPath(params.file_ref_gzip)
     if(fileexist_b("${params.file_ref_gzip}.csi")){
     rs_infogz=file_ref_gzip.combine(channel.fromPath(params.file_ref_gzip+'.csi',checkIfExists:true))
   }else {
     bcftools_index_vcf(file_ref_gzip)
     rs_infogz=bcftools_index_vcf.out
   }
   pp['rs_info']  = rs_infogz
 }
 return(pp) 

}

workflow format_plink_invcf {
 take:
   plink
   outputdir
   outputpat
   fasta
  main :
   pp=check_params(plink, fasta)
   if(!(pp['rs_info']!=null)) {
   updateplk_rsname_norsinfo(pp['bfile'],  channel.of("$outputdir/rs_update"))
   }else{
     updateplk_rsname(pp['rs_info'], pp['bfile'], channel.of("$outputdir/rs_update"))
     deletedmultianddel(updateplk_rsname.out.plink)
     refallele(deletedmultianddel.out, updateplk_rsname.out.rs_info)
     if(params.convertinvcf_parralchro==0){
      //    tuple path(bed), path(bim), path(fam), path(gz_info), path(fast), path(fastaindex),val(chro)
      convertInVcf(refallele.out.combine(getparams.out.rs_infogz).combine(getparams.out.fasta).combine(channel.of("")).combine(channel.of("$outputdir/")))
      checkfixref(convertInVcf.out.vcf, getparams.out.fasta, channel.of("$outputdir/fixref"))
      checkVCF(convertInVcf.out.vcf, getparams.out.fasta, channel.of("$outputdir/checkVCF"))
      vcf = checkVCF.out
     }else{
      list_chro(refallele.out)
      list_chro2=list_chro.out.flatMap { list_str -> list_str.split()}
      convertInVcf(refallele.out.combine(getparams.out.rs_infogz).combine(getparams.out.fasta).combine(list_chro2).combine(channel.of("$outputdir/vcf/")))
      mergevcf(convertInVcf.out.vcf.collect(), channel.of("$outputdir/"))
      vcf= mergevcf.out
    }
   }
   if(pp[['rs_info']]!=null)get_sex_forsanger(getparams.out.plink, vcf, "$outputdir/") 
  emit : 
      vcf= vcf
}

