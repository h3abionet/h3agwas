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

//Process 1: imported the plink files (bed, bim, fam) from raw_src channel and
//exported them through 4 new channels: out_ch, ch_select_rs_format, 
//fastlmm_assoc_ch, rel_ch_fastlmm

process aka{
    echo true

    input: 
      set file(bed), file(bim), file(fam) from raw_ch

    output:
      set bed, bim, fam into (out_ch, ch_select_rs_format, fastlmm_assoc_ch, rel_ch_fastlmm, fam_ch_fast)

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

num_assoc_cores = params.mperm == 0 ? 1 : Math.min(10,params.max_plink_cores)

supported_tests = ["assoc","fisher","model","cmh","linear","logistic"]

requested_tests = supported_tests.findAll { entry -> params.get(entry) }

covariate = ""
gotcovar  = 0
pheno     = ""

balise_filers_rel=1

//Relatedness

//Process 4: Checking if any of the following parameters are provided: boltlmm, gemma, 
//fastlmm or fastgwa (1 or 0). If any one is provided, this condition will be met. 
//In condition 2, 

process select_rs_format {

    cpus max_plink_cores
    memory plink_mem_req
    time   params.big_time
    echo true

    input:
       set file(bedd),file(bimm), file(famm) from ch_select_rs_format

    output:
       set file("${params.sampling}.prune.in") into  filers_matrel_mat_fast, filers_matrel_mat_GWA, filers_matrel_mat_gem, filers_matrel_bolt, filers_count_line

    script:
        base = params.input_testing
        prune = "${params.output_testings}-prune"

        """
        plink --bfile ${base} --indep-pairwise ${params.sample_snps_rel_paramplkl} --out $prune   --threads ${params.max_plink_cores}
        plink --bfile ${base} --indep-pairwise ${params.sample_snps_rel_paramplkl} --out ${params.sampling}   --threads ${params.max_plink_cores}
        """
}

balise_filers_rel=0

filers_matrel_mat_fast=file('NO_FILE')
filers_matrel_mat_GWA=file('NO_FILE')
filers_matrel_mat_gem=file('NO_FILE')

process buildBoltFileSnpRel{

    memory params.bolt_mem_req
    time   params.big_time

    input:
        set file(bed),file(plinksbim), file(fam) from src_ch

    output :
        file(output) into filers_matrel_bolt_else

    script :
        output=plinksbim.baseName+".rs.choice"
        start = "start"
        
        """
        shuf -n 950000 $plinksbim | awk '{print \$2}' > $output
        shuf -n 950000 $plinksbim | awk '{print \$2}' > $start
        echo 7
        """
}

// filers_matrel_mat_fast=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
// filers_matrel_bolt=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
// filers_matrel_mat_GWA=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)
// filers_matrel_mat_gem=Channel.fromPath(params.file_rs_buildrelat, checkIfExists:true)

if (params.data != "") {

  params.data = file(params.data)

  if (params.covariates != "") {
      gotcovar = 1
  }
}

  process extractPheno {

    input:
     set file(data) from data_ch

    output:
     file("pheno.phe") into pheno_ch

    script:
    all_phenos = params.covariates.length()>0 ? params.pheno+","+params.covariates : params.pheno
      """
      extractPheno.py $data ${all_phenos} "pheno.phe"
      """
}

pheno_label_ch = Channel.from(params.pheno.split(","))

process showPhenoDistrib {

    input:
    file(data) from data_ch

    output:
      file ("B050*") into pheno_report_ch

    script:
      "phe_distrib.py --pheno ${params.pheno} $data B050 "
}

data_ch_fastlmm = Channel.fromPath(params.data, checkIfExists:true)

if (params.covariates)
     covariate_option = "--cov_list ${params.covariates}"
  else
     covariate_option = ""

// fam_ch_fast = Channel.from(params.input_testing_fam)

  process  getFastLmmPhenosCovar {

    input:
      file(covariates) from data_ch_fastlmm
      set file(fam) from fam_ch_fast

    // output:
    //   set file(phef), file(covfile) into fastlmm_data_ch
    //   stdout into pheno_cols_ch_fastlmm

    script:
      
      """ echo 9"""
    //   base = fam.baseName
      phef = "sampleA_fastlmm_n.phe"
      covfile = "sampleA_fastlmm_n.cov"

      """
       all_covariate.py --data  $covariates --inp_fam  $fam $covariate_option \
             --pheno ${params.pheno} --phe_out ${phef}  --cov_out $covfile --form_out 3
       """

  }