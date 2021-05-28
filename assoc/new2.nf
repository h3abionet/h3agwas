#!/usr/bin/env nextflow
nextflow.enable.dsl=2

plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores 

src_ch = Channel.from(params.input_testing)

src_cha = Channel.fromPath(params.input_testinga)
src_chb = Channel.fromPath(params.input_testingb)
src_chc = Channel.fromPath(params.input_testingc)

workflow {

    out_ch = importplink(src_cha.combine(src_chb).combine(src_chc))

    pca_out_ch = computePCA(out_ch)

    report_pca_ch = drawPCA(pca_out_ch)
}


//Process 1: imported the plink files (bed, bim, fam) from raw_src channel and
//exported them through 4 new channels: out_ch, ch_select_rs_format, 
//fastlmm_assoc_ch, rel_ch_fastlmm

process importplink{
    echo true

    input: 
        tuple file(bed), file(bim), file(fam)

    output:
        tuple file(bed), file(bim), file(fam)

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

    cpus max_plink_cores
    memory plink_mem_req
    time   params.big_time

        input:
        tuple file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam')

        output:
        tuple file("${params.sampled}.eigenval"), file("${params.sampled}.eigenvec")

        script:

        base = "cleaned"
        prune = "${params.output_testing}-prune"

        """
        plink --bfile ${params.input_testing} --indep-pairwise 100 20 0.2 --out check
        plink --keep-allele-order --bfile ${params.input_testing} --extract check.prune.in --make-bed --out $prune
        plink --threads $max_plink_cores --bfile $prune --pca --out ${params.output_testing} 
        plink --threads $max_plink_cores --bfile $prune --pca --out ${params.sampled}
        """
}

//Process 3: Imported the eigen files from the pca_out_ch and used the drawPCA.py script to
//to plot the PCA graph.

process drawPCA {
    input:
        tuple file(eigvals), file(eigvecs)
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