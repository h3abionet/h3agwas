#!/usr/bin/env nextflow
nextflow.enable.dsl=2

src_cha = Channel.fromPath(params.input_testinga)
src_chb = Channel.fromPath(params.input_testingb)
src_chc = Channel.fromPath(params.input_testingc)

bed = Channel.fromPath("${params.input_testing}.bed")
bim = Channel.fromPath("${params.input_testing}.bim")
fam = Channel.fromPath("${params.input_testing}.fam")

// bed = Paths.get("${params.input_testing}.bed").toString()
// bim = Paths.get("${params.input_testing}.bim").toString()
// fam = Paths.get("${params.input_testing}.fam").toString()

// process importplink{
//     echo true

//     input: 
//         tuple file(bed), file(bim), file(fam) from bed,bim,fam

//     // output:
//     //     tuple bed, bim, fam into (out_ch)

//     script:
        
//         // bed = file(bed)
//         // bim = file(bim)
//         // fam = file(fam)

//         """
//         echo transformed files correctly!
//         """
// }

process foo {
input:
  tuple X, 'some-file.bam'
 script:
   '''
   echo 99999
   '''
}