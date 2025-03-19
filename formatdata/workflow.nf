include {bcftools_index_vcf} from '../modules/process_vcf.nf'
include {list_chro} from '../modules/utils_plink.nf'
include {updateplk_rsname;updateplk_rsname_norsinfo;deletedmultianddel;refallele;convertInVcf;checkfixref;checkVCF;mergevcf;get_sex_forsanger} from './process.nf'
include {fileexist_b;fileexist_param;fileexist} from '../modules/fct_groovy.nf'
include {dl_fasta_wf} from '../modules/dl_wf.nf'

workflow getparams{
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

workflow format_plink_invcf {
 take:
   plink
   outputdir
   outputpat
   fasta
  main :
  //spaces/jeantristan/agora/imputeddata/all/job-20241209-124900-327/local/chr11.dose.vcf.gz, /spaces/jeantristan/agora/imputeddata/convert_and_merge/all/work/ba/d6e6622f6774531f1f806096afa73a/hg38.fa_clean.fa.gz, 0.01, 0.008, 0, -1, .//convertvcf/cleanvcf]
   getparams(plink, fasta, vcf_imputeformat)
   if(!(getparams.out.rs_infogz!=null)) {
   updateplk_rsname_norsinfo(getparams.out.plink,  channel.of("$outputdir/rs_update"))
   }else{
     updateplk_rsname(getparams.out.rs_infogz, getparams.out.plink, channel.of("$outputdir/rs_update"))
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
   if(params.sexinfo_available)get_sex_forsanger(getparams.out.plink, vcf, "$outputdir/") 
  emit : 
      vcf= vcf
}

