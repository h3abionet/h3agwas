#!/usr/bin/env nextflow

import java.nio.file.Paths

ch_select_rs_format=Channel.create()

data_ch = file(params.data)

params.max_plink_cores = 4
params.plink_mem_req = '6GB'
params.other_process_mem_req = '10G'
params.big_time             = '1000h'

plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores 

src_ch = Channel.from(params.input_testing)
raw_ch = Channel.from(params.input_testing)

rel_ch_fastlmm = Channel.create()
fastlmm_assoc_ch = Channel.create()

bed = Paths.get("${params.input_testing}.bed").toString()
bim = Paths.get("${params.input_testing}.bim").toString()
fam = Paths.get("${params.input_testing}.fam").toString()

// Process 1: imported the plink files (bed, bim, fam) from raw_src channel and
// exported them through 4 new channels: out_ch, ch_select_rs_format, 
//fastlmm_assoc_ch, rel_ch_fastlmm

process aka{
  echo true

  input: 
    set file(bed), file(bim), file(fam) from raw_ch

  output:
    set bed, bim, fam into (out_ch, ch_select_rs_format, fastlmm_assoc_ch, rel_ch_fastlmm)

  script: 
    bed = file(bed)
    bim = file(bim)
    fam = file(fam)

    """
    echo transformed files correctly!
    """
}

// Process 2: Computing PCA using the imported plink files from the out_ch channel 
//and exported two files, the eigenvalues and the eigenvector files 
//into the pca_out_ch channel.

process computePCA {
  echo true
  publishDir params.output, overwrite:true, mode:'copy'
  cpus max_plink_cores
  memory plink_mem_req
  time   params.big_time

  input:
    set file('cleaned.bed'),file('cleaned.bim'),file('cleaned.fam') from out_ch

  output:
    set file("${params.sampled}.eigenval"), file("${params.sampled}.eigenvec") into pca_out_ch

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
    set file(eigvals), file(eigvecs) from pca_out_ch
  output:
    set file (output), file ("B040-pca.tex") into report_pca_ch
    publishDir params.output, overwrite:true, mode:'copy',pattern: "*.pdf"
  script:
    base=eigvals.baseName
    cc_fname = 0
    cc       = 0
    col      = 0

    output="${params.first}-pca.pdf"
    template "drawPCA.py"
}