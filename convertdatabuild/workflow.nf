include {dl_fasta_wf} from '../modules/dl_wf.nf'                                               
include {picard_liftover} from './process.nf'

workflow getparams {
  take :
    vcf_qc
    ref_in
    ref_out
  main :
   if(params.vcf=='' && params.vcf_list==''){
       if(vcf_qc==null){
         println "no vcf :"
         System.exit(1)
       }
  }else{
   if(params.vcf!=''){
     vcf_qc = channel.fromPath(params.vcf)
  }else {
    vcf_qc=Channel.fromPath(file(params.vcf_list).readLines(), checkIfExists:true)
  }
 } 
 dl_fasta_wf(params.fasta_convert,ref_out, params.ftp_fasta_convert, 'picard')
 if(params.crossmap_chains!=''){
   crossmap=channel.fromPath(params.crossmap_chains, checkIfExists:true)
 }
 emit :
   vcf = vcf_qc
   fasta = dl_fasta_wf.out.fasta_index
   fasta_picard = dl_fasta_wf.out.fasta_picard
   crossmap = crossmap
}

workflow crossmap_vcf{
  take :
     vcf_qc
  main : 
    getparams(vcf_qc, params.build_genome, params.build_genome_convert)
    out=getparams.out.vcf.flatMap{it->it.Name.toString().replace('.vcf.gz', '.'+params.build_genome_convert+'.vcf.gz')}
    outreject=getparams.out.vcf.flatMap{it->it.Name.toString().replace('.vcf.gz', 'reject.vcf.gz')}
    picard_liftover(getparams.out.vcf.combine(getparams.out.crossmap).combine(getparams.out.fasta_picard).combine(out).combine(outreject).combine(channel.of("$params.output_dir/convert_build${params.build_genome_convert}")))
    emit :
       vcf = picard_liftover.out.vcf
}
