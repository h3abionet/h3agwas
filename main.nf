#!/usr/bin/env nextflow                                                         
                                                                                
/*                                                                              
 * Authors       :                                                              
 *                                                                              
 *      Jean-Tristan Brandenburg                                                
 *      Scott Hazelhurst                                                        
 *                                                                              
 *  On behalf of the H3ABionet Consortium                                       
 *  2015-2024
 *                                                                              
 *                                                                              
 * Description  : Nextflow pipeline for association .                              
 *                                                                              
 *(C) University of the Witwatersrand, Johannesburg, 2016-2024
 *    on behalf of the H3ABioNet Consortium                                     
 *This is licensed under the MIT Licence. See the "LICENSE" file for details    
 */


# quality control 
include { qc } from "./qc/workflow.nf"
include { checkresume } from "./modules/fct_groovy.nf"
//include { assoc} from "./assoc/assoc.nf"
//nextflow.enable.moduleBinaries = true

workflow.onComplete = {
    // any workflow property can be used here
    println "Pipeline complete"
    println "Command line: $workflow.commandLine"
}

workflow {
  checkresume()
  //if( !nextflow.version.matches('21.04+') ) {
  //  println "This workflow requires Nextflow version 21.04 or greater -- You are running version $nextflow.version"
  //  exit 1
  //}
  if (params.qc == 1 || params.qc) {
        qc()
  }
}
