plink_mem_req = params.plink_mem_req
other_mem_req = params.other_process_mem_req
max_plink_cores = params.max_plink_cores 

process importplink{
    echo true
    tag "Importing plink files"

    input: 
        path bed
        path bim
        path fam

    output:
        path bed
        path bim
        path fam

    script:
        bed = bed
        bim = bim
        fam = fam

        """
        echo transformed files correctly!
        """
}

process computePCA {
    tag "Converting to eigens"
    publishDir "${params.output}/Eigens", overwrite:true, mode:'copy'

    cpus max_plink_cores
    memory plink_mem_req
    time   params.big_time

        input:
            path bed
            path bim
            path fam

        output:

            path "${params.sampledd}.eigenval"
            path "${params.sampledd}.eigenvec"

        script:
            prune = "${params.output_testers}-prune"

            """
            plink --bfile ${params.input_testers} --indep-pairwise 100 20 0.2 --out check
            plink --keep-allele-order --bfile ${params.input_testers} --extract check.prune.in --make-bed --out $prune
            plink --threads $max_plink_cores --bfile $prune --pca --out ${params.output_testers} 
            plink --threads $max_plink_cores --bfile $prune --pca --out ${params.sampledd} 
            """
}

process drawPCA {
    tag "Drawing PCA graphs"
    publishDir "${params.output}/PCA", overwrite:true, mode:'copy'

    input:
        path eigvals
        path eigvecs

    output:
        path "sampleA-pca.pdf"
        path "eigenvalue.pdf"
        path "B040-pca.tex"

    script:
        base=eigvals.baseName
        cc_fname = 0
        cc       = 0
        col      = 0

        output="sampleA-pca.pdf"
        output_2="eigenvalue.pdf"
        output_3="B040-pca.tex"
        template "drawPCA.py"
}
