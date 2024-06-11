include {bcftools_index_vcf} from '../modules/vcf.nf'
include {list_chro} from '../modules/utils_plink.nf'
include {updateplk_rsname;updateplk_rsname_norsinfo;deletedmultianddel;refallele;convertInVcf;checkfixref;checkVCF;mergevcf} from './process.nf'
include {fileexist_b;fileexist_param;fileexist} from '../modules/fct_groovy.nf'
include {df_fasta_wf} from '../modules/dl_wf.nf'

workflow getparams {
  take :
   plink
  main :
   if(plink==null){
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
 df_fasta_wf(params.fasta, params.build_genome, params.ftp_fasta)
 emit :
  rs_infogz= rs_infogz
  fasta = df_fasta_wf.out.fasta_index
  plink = plink
}

workflow format_plink_invcf {
 take:
   plink
   outputdir
   outputpat
  main :
   getparams(plink)
   if(getparams.out.rs_infogz==null) {
   updateplk_rsname_norsinfo(getparams.out.plink,  channel.of("$outputdir/rs_update"))
   }else{
     updateplk_rsname(getparams.out.rs_infogz, getparams.out.plink, channel.of("$outputdir/rs_update"))
     deletedmultianddel(updateplk_rsname.out.plink)
     refallele(getparams.out.plink, updateplk_rsname.out.rs_info)
     if(params.convertinvcf_parralchro==0){
      //    tuple path(bed), path(bim), path(fam), path(gz_info), path(fast), path(fastaindex),val(chro)
      convertInVcf(refallele.out.combine(getparams.out.rs_infogz).combine(getparams.out.fasta).combine(channel.of("")).combine(channel.of("$outputdir/")))
      checkfixref(convertInVcf.out.vcf, getparams.out.fasta, channel.of("$outputdir/fixref"))
      checkVCF(convertInVcf.out.vcf, getparams.out.fasta, channel.of("$outputdir/checkVCF"))
      vcf = checkVCF.out
     }else{
      list_chro(getparams.out.plink)
      list_chro2=list_chro.out.flatMap { list_str -> list_str.split()}
      convertInVcf(getparams.out.plink.combine(getparams.out.rs_infogz).combine(getparams.out.fasta).combine(list_chro2).combine(channel.of("$outputdir/vcf/")))
      mergevcf(convertInVcf.out.vcf.collect(), channel.of("$outputdir/"))
      vcf= mergevcf.out
    }
   }
  emit : 
      vcf= vcf
}
