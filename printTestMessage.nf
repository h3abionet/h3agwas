/* 
* PRINT TEST MESSAGE 
* =================
*/

/************************************************/

nextflow.enable.dsl=2

include {
    getChannels;
    getInputChannels;
    printToScreen;
} from "${projectDir}/modules/testMessage.nf"

workflow {

    message = getChannels()
    messages = getInputChannels()

    printToScreen(message, messages) | view()

}
