include {join2channel;check_pheno_bin} from './all.nf'
include { strmem } from "../../modules/fct_groovy.nf"                           

                                                                                
                                                                                
process computeTest {                                                        
      // template                                                               
      memory { strmem(params.high_memory) + 5.GB * (task.attempt -1) }            
      cpus params.max_cpus                                                        
      maxRetries 10                                                               
      time params.big_time                                                      
      input:                                                                    
       tuple val(test), val(pheno),val(pheno_type),val(covariates), path(bed), path(bim), path(fam), path(phenof) 
      publishDir "${params.output_dir}/assoc/${test}", overwrite:true, mode:'copy'    
      output:                                                                   
        tuple val(test), val(pheno), path("${outfname}.assoc.*"),emit : plink
        tuple path('*.log'),path("*.report"),emit : report 
      script:                                                                   
       println(test+' '+pheno_type)
       base = bed.baseName
       perm = (params.mperm == 0 ? "" : "--mperm=${params.mperm}")              
       adjust = (params.adjust ? "--adjust" : "")                               
       outfname = "${pheno}"                                               
       pheno_cmd = "--pheno $phenof --pheno-name $pheno"              
       covarite=""
       if (covariates!='') covariate = "--covar ${phenof} --covar-name ${covariates} "
       sexallow="--allow-no-sex"                                               
       """
       plink --bfile $base $covariate $pheno_cmd --threads ${params.max_cpus} --${test} $perm $adjust --out $outfname $sexallow
      cp .command.sh "${pheno}"_"$test".cmd.report                          
      cp .command.log "${pheno}"_"$test".log.report                         
      cp .command.err "${pheno}"_"$test".err.report  
       """
}   
process  format_sumstat_plink{
  input :
  

}
workflow plink {
  take :                                                                        
   data                                                                         
   pheno                                                                        
   pheno_bin                                                                    
   covariates                                                                   
   covariates_type                                                              
   plink                                                                        
   plink_rel                                                                    
   vcf                                                                          
   vcf_balise                                                                   
   bgen                                                                         
   bgenlist                                                                     
   bgen_sample                                                                  
   bgen_balise                                                                  
   bgenlist_balise                                                              
   listchro                 
 main :                                                                         
   check_pheno_bin(pheno,pheno_bin,data,'plink')
   supported_tests = ["assoc","fisher","model","cmh","linear","logistic"]          
   requested_tests = channel.from(supported_tests.findAll { entry -> params.get(entry) }       )
   println(requested_tests)
   phenol = pheno.flatMap { list -> list.split(',') }  // Groovy-style lambda for splitting
   pheno_binl = pheno_bin.map { list -> list.split(',') }.flatMap { it.size() == 1 ? it.collect { it } * npheno : it }
   join2channel(phenol,pheno_binl) 
   alltest=requested_tests.combine(join2channel.out).combine(covariates).combine(plink).combine(check_pheno_bin.out.data)
   alltest=alltest.filter{t-> ((t[2]=='1' && (t[0] in ['fisher','logistic'])) || (t[2]=='0' && (t[0] in ['linear','assoc'])))}
   computeTest(alltest) 
   //format_plink(computeTest.out )
}
