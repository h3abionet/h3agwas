include {MD5_plk as MD5_in} from '../modules/utils.nf'                          
include {is_nullfile; fileexist} from '../modules/fct_groovy.nf'                

workflow check_params{                                                          
take :                                                                         
  data
  bfile
 main :                                                                        
  filescript=file(workflow.scriptFile)                                         
  projectdir="${filescript.getParent()}" 

  if (workflow.repository)                                                      
    wflowversion="${workflow.repository} --- ${workflow.revision} [${workflow.commitId}]"
  else                                                                          
    wflowversion="A local copy of the workflow was used"   

 if(data){                                                                      
  phenotype_ch = data                                                           
 }else {     


 }

}


workflow assoc {
 take :
   data
   bfile
   vcf
   bimbam
   bgen
  println "In ASSOC top"
 // add pcs

       println "In ASSOC middle"
}
