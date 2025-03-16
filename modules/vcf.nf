include {list_chro;splitvcf;splitvcf2} from './process_vcf.nf'
include {strmem} from './fct_groovy.nf'
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


process buildindex {
   memory { strmem(params.low_memory) + 5.GB * (task.attempt -1) }
   errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
   maxRetries 3
    label 'utils'
    input :
       path(vcf) 
    output :
        tuple path(vcf), path(vcfindex) 
    script :
       vcfindex=vcf+".csi"
       """
       tabix -C -p vcf $vcf
       """
}

process getlistchro{
 label 'utils'
 input :
    path(vcf)
 output :
   stdout
 script :
   """
   zcat $vcf|grep -v "#"|awk '{print \$1}' |uniq|sort|uniq
   """
}
