#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {

    importplink; 
    computePCA; 
    drawPCA

} from './modules/PCA.nf';

workflow {

    src_cha = Channel.fromPath(params.first)
    src_chb = Channel.fromPath(params.second)
    src_chc = Channel.fromPath(params.third)

    importplink(src_cha, src_chb, src_chc)

    computePCA(importplink.out)

    drawPCA(computePCA.out)

}