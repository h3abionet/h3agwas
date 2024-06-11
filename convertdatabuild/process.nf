include {strmem}  from '../modules/fct_groovy.nf'

process picard_fastaindex{
 memory { strmem(params.high_memory) * task.attempt }                           
 errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }          
 maxRetries 3                                                                   
 label 'picard'                                                                 
 input :                                                                        
   path(fasta)
   val(outputdir)
publishDir "${outputdir}/",  mode:'copy'                                        
 output :                                                                       
   tuple path("$fasta"), path("$out")
 script :                                                                       
   out=fasta+'.dict'
   """                                                                          
   java -jar ${params.bin_picard} CreateSequenceDictionary  R=$fasta  O=$out
   """   
}

process  picard_liftover {
 memory { strmem(params.high_memory) * task.attempt }                                   
 errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }         
 maxRetries 3   
 label 'picard'
 input :
   tuple path(vcf),path(chain), path(fasta),path(fasta2),val(vcfout), val(outputreject), val(outputdir)
publishDir "${outputdir}/",  mode:'copy'
 output :
   path("$vcfout"), emit: vcf
   path("$outputreject"), emit: vcf_error
 script :
   """
   java -jar ${params.bin_picard} LiftoverVcf  -I $vcf -O $vcfout -CHAIN $chain -REJECT $outputreject -R $fasta --MAX_RECORDS_IN_RAM ${params.picard_max_record}
   """
}
