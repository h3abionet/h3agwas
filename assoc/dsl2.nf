#!/usr/bin/env nextflow
nextflow.enable.dsl=2

src_cha = Channel.fromPath(params.a)
src_chb = Channel.fromPath(params.b)
src_chc = Channel.fromPath(params.c)

// src_cha.combine(src_chb).combine(src_chc)

//Process 1: imported the plink files (bed, bim, fam) from raw_src channel and
//exported them through 4 new channels: out_ch, ch_select_rs_format, 
//fastlmm_assoc_ch, rel_ch_fastlmm

process importplink{
    echo true

    input: 
        set file(bed), file(bim), file(fam)

    output:
        set  bed, bim, fam

    script:
        
        bed = file(bed)
        bim = file(bim)
        fam = file(fam)

        """
        echo transformed files correctly!
        """
}

//Process 2: Computing PCA using the imported plink files from the out_ch channel 
//and exported two files, the eigenvalues and the eigenvector files 
//into the pca_out_ch channel.

process computePCA {
    echo true
    publishDir params.output, overwrite:true, mode:'copy'

    cpus params.max_plink_cores 
    memory params.plink_mem_req
    time   params.big_time

        input:
        tuple file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam')
        output:
        tuple file("${params.sampled}.eigenval"), file("${params.sampled}.eigenvec")

        script:

        """
        echo believe!
        """

        prune = "${params.output_testing}-prune"

        """
        plink --bfile ${params.input_testing} --indep-pairwise 100 20 0.2 --out check
        plink --keep-allele-order --bfile ${params.input_testing} --extract check.prune.in --make-bed --out $prune
        plink --threads $max_plink_cores --bfile $prune --pca --out ${params.output_testing} 
        plink --threads $max_plink_cores --bfile $prune --pca --out ${params.sampled}
        """
}

// //Process 3: Imported the eigen files from the pca_out_ch and used the drawPCA.py script to
// //to plot the PCA graph.

process drawPCA {
    input:
        tuple file(eigvals), file(eigvecs) from pca_out_ch
    output:
        tuple file (output), file ("B040-pca.tex")
        publishDir params.output, overwrite:true, mode:'copy',pattern: "*.pdf"
    script:
        base=eigvals.baseName
        cc_fname = 0
        cc       = 0
        col      = 0

        output="${params.first}-pca.pdf"
        template "drawPCA.py"
}

workflow {

  importplink(src_cha.combine(src_chb).combine(src_chc))
  computePCA(importplink.out)
  drawPCA(computePCA.out)

}