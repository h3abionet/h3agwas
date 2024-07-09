include {list_chro;splitvcf} from './process_vcf.nf'
workflow getparams {
 take :
  vcf 
  main :
  if(!vcf){
    vcf=channel.fromPath(params.vcf)
  }
 emit : 
   vcf = vcf
}

workflow split_vcf{
  take :
   vcf
   outputpat
   outputdir
  main :
   getparams(vcf)
   list_chro(getparams.out.vcf) 
   splitvcf(getparams.out.vcf, channel.of(outputdir),channel.of(outputpat))
  emit :
   vcf=splitvcf.out
}
