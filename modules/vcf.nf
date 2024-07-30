include {list_chro;splitvcf;splitvcf2} from './process_vcf.nf'
workflow getparams {
 take :
  vcf 
  main :
  if(!vcf){
    vcf=channel.fromPath(params.vcf, checkIfExists:true)
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
   list_chro2=list_chro.out.flatMap { list_str -> list_str.split()}
   splitvcf2(getparams.out.vcf.combine(list_chro2).combine(channel.of(outputdir)).combine(channel.of(outputpat)))
  emit :
   vcf=splitvcf2.out
}
