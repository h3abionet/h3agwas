<img src="../auxfiles/H3ABioNetlogo2.jpg"/>
#Heritabilies estimation  : `heritabilities/main.nf`

This section describes a pipeline in devlopment, objectives is estimated heritabilities with various way, we developped : ldlc, grmel of bolt and greml of gcta, gemma
two distincs approaches should be considered :
  * based on relatdness matrix and phenotype as gcta, bolt, gemma
  * based on gwas result as implemented in ldlc and gemma


### Installation
need python3, gcta, ldlc, bolt and gemma, R with ggplot2

### Running
The pipeline is run: `nextflow run assoc/esth2-assoc.nf`

The key options are:
* `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
* `output_dir` : output directory
* `output_pat` : output pattern
* ̀ data` : same option that _assoc/plink-assoc.nf_, file is optional, used for gemma, bolt and gcta
  * `pheno` :phenotypes used in data to computed in gemma, bolt
* `file_gwas` : one ore more one file gwas, used for ldsc and gemma, to defined column of files :
  * ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
  * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
  * `head_rs` : rs header column [default : "SNP"]
  * `head_beta` : beta header colum [default : "BETA"]
  * `head_se`  : column for standard error of beta "SE"
  * `head_A1` : column for A0 :[default : "ALLELE0" ]
  * `head_A2` : column for A0 :[default : "ALLELE2" ]
  * `head_freq` : freq header [ default : A1Freq],
  * `head_n`: N header, used just for ldsc, if not present, `Nind` must be initialize.
  * `Nind` : if `head_n` not initialise, must be initialise, individuals number for each gwas file, separate by comma
* `ldsc_h2` : need a estimation of h2 by ldc : 1 [default : 0]:
  * [LDSC](https://github.com/bulik/ldsc) computes heritabilies between gwas value using LD information.
  * `ldsc_bin` : binary for ldsc
  * `dir_ref_ld_chr` : folder containing ld information, 1 file by begin by chromosome num without chr, see : `--ref-ld-chr` in ldsc manual
  * ̀`ldsc_mem_req`
  * `ldsc_h2_multi` : computing genetic correlation. between different gwas result
  * `munge_sumstats_bin` : binary for munge
  * `ldsc_h2opt` : other option for ldsc
  * output :
* `gemma_h2` : estimation using gemma,
  * `gemma_bin` : gemma binary [ default "gemma"]
  * `gemma_h2` : do you want a estimation with gemma of heritabilies with relatdness matrix [default : 0] :
  * `gemma_h2_pval` : do you want a estimation of heritabilities with gemma using p-value [ default : 0]
   * need file of pvalue see `file_gwas`
   * Z obtained with beta/se
  * `gemma_h2_typeest` do you wang a (1) or reml (2)[default : "1"]
  * `gemma_mat_rel` for `gemma_h2`, have you a pre estimated relatdness otherwise it wil be computed with plink file [default : none]
  * `gemma_num_cores` : Core number [ default : 6]
  * `gemma_mem_req` : memory for gemma [ default : "10GB" ]
  * `gemma_relopt` = 1
  * file_rs_buildrelat : list of rs used to build relatdness [default : None ]
* `bolt_h2` : estimation using bolt
  *  see [manual](https://data.broadinstitute.org/alkesgroup/BOLT-LMM/)
  * file_rs_buildrelat : list of rs used to build relatdness [default : None ]
     * if SNPs is higher than 950,000, 950,000 SNPs are chosen randomly to build the model (see --modelSnps option in bolt)
  * `bolt_covariates_type` : for bolt need to define if covariate is binary (0) or continue (1), a comma-separated list as same order than covariates
  * `bolt_ld_score_file` : A table of reference LD scores for boltlmm is needed to calibrate the BOLT-LMM statistic (option in boltlmm --LDscoresFile),to choice a specific column in Ld file you can use `bolt_ld_scores_col` option (by default : LDSCORE) if option is not provided --LDscoresUseChip used.
  * `bolt_num_cores` if bolt is used set this up to 8
  * `bolt_mem_req` memory required for boltlmm, (default : 6GB)
  * `bolt_h2_multi`  : to have multi variance between traits [ default : 0]
* `gcta` :
  * general option :
   * `output_dir` : plink directory
   * `output_pat` : basename plink
   * `data` : data phenotypes
   * `pheno` : pheno to estimate h2
   * `covariates` : covariate to estimate heritabilities
   * `gcta_mem_req` : [default 20GB]
   * `gcta_bin` : binary for gcta
   * `gcta_h2_ldscore` [default 200kb]
   * `gcta_h2_multi` : computed co heritability between phenotype
   * `gcta_h2_imp` : 0 : version for low density of snps as dna chip
   * `gcta_h2_imp` : 1 version for high density of snps as data imputed, described [here](https://cnsgenomics.com/software/gcta/#GREMLinWGSorimputeddata) with au
   * `gcta_h2_imp` : 1 version for high density of snps as data imputed, described [here](https://cnsgenomics.com/software/gcta/#GREMLinWGSorimputeddata) with au
    * `gcta_h2_grmfile` : extension of grm output for 3 files `gcta_h2_grmfile`.grm.bin, `gcta_h2_grmfile`.grm.id `gcta_h2_grmfile`.grm.N.bin otbained with command line : gcta64 --bfile plk --make-grm --out outhead, if `gcta_h2_grmfile` is empty default, pipeline will generate files and output will be in gctagrm folder, __if not present will be generated__ :
    * for imputed data, we used algorithm of [here](https://cnsgenomics.com/software/gcta/#GREMLinWGSorimputeddata) and build :
     * `grm_cutoff` : `grm-cutoff` in gcta (alias of `--rel-cutoff`) If used in conjunction with a later calculation (see the order of operations page for details), --rel-cutoff excludes one member of each pair of samples with observed genomic relatedness greater than the given cutoff value (default 0.025) from the analysis. Alternatively, you can invoke this on its own to write a pruned list of sample IDs to plink.rel.id. [default : 0.025]
      * Step 1: segment based LD score
      * Step 2 : stratify the SNPs by LD scores of individual SNPs in R
      * Step 3: making GRMs using SNPs stratified into different groups
      * steps was very long you can provide `gcta_h2_grmfile` option to skip with already snps stratified
      * `gcta_mem_reqmgrm` : option memory for the step
    * `gcta_h2_mgrmfile` : if file not provide pipeline will do step described before[default None]
   * `gcta_mem_reqmgrm` : [default 40GB]
   * `params.gcta_reml_alg` : see reml-alg :  Specify the algorithm to run REML iterations, 0 for average information (AI), 1 for Fisher-scoring and 2 for EM. The default option is 0, i.e. AI-REML, if this option is not specified.  [oa]

