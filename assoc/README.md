<img src="../auxfiles/H3ABioNetlogo2.jpg"/>

#  Association pipeline: `assoc/main.nf`

This workflow has been extensively expanded by Jean-Tristan Brandenburg

An association study is a complex analysis and each analysis has to consider
* the disease/phenotype being studied and its mode of inheritance
* population structure
* other covariates

For this reason it is difficult to build a high quality, generic pipeline to do an association study. 

The purpose of this pipeline is to perform a very superficial initial analysis that can be used as one piece of information to guide a rigorous analysis. Of course, we would encourage users to build their own Nextflow script for their rigorous analysis, perhaps using our script as a start.

Our script, *assoc* takes as input PLINK files that have been through quality control and 
* does a principal component analysis on the data, and produces pictures from that; 
* performs a simple association test giving odds ratio and  raw and adjusted _p_ values

## 1. Running the main association testing

The pipeline is run: `nextflow run assoc`

The key options are:
* `input_dir`, `output_dir`: where input and output goes to and comes from;
* `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
* `data`: a tab-separated file that contains phenotype data for your particpants. The row is a header line, with one participant per line after that. The first two columns should FID IID followed by any phenotype values you want to use for testing or for covariates. It can contain other data too -- as long as the ones that you need are in this file.
* `pheno`: a comma-separated list of phenotypes that you want to test. Each phenotype is tested separately. If you wish to do a transform to the phenottype, you can suffix the phenotype with _/function_ where function is a numpy function. For example you might say `--pheno bmi/np.log` which will apply the log function to the phenotype. Any numpy function can be used, typical uses being np.log and np.sqrt. We plan to support user provision of a user-given function.
* `covariates`: a comma-separated list of phenotypes that you want to use
* `exclude_snps` option to exclude some snps active just for boltlmm (see `--exclude` in boltlmm manual) : SNP ID must be first column (default none)
*  `print_pca` : by default pipeline compute and print pca (`print_pca`=1), if you want avoid this step (`print_pca` = 0)
*  `file_rs_buildrelat` : file with rs list (one by lines) to build genetics models (relatdness), for gemma `-snps` for boltlmm `--modelSnps`
* `genetic_map_file` : genetic map used in boltlmm 

By default a chi2 test for association is done. But you can do multiple different tests in one run by settintg the appropriate parameter to 1. Note at least one must be set to 1

 * `assoc` : should a chi2 test be used (0 or 1)
 * `fisher`: Fisher exact test (default 0)
 *  `linear`: should linear regreession be used?  (default 0)
 *  `logistic`: should linear regression be used? (default 0)
 *  `gemma`: should gemma be used? (default 0)
    *  see [manual](www.xzlab.org/software/GEMMAmanual.pdf)
    *  `gemma_num_cores`: if gemma is used set this up to 8
    *  `gemma_mem_req`: For 10k samples, 2 million SNPs, we needed 4GB of RAM (default : 6GB)
    *  `gemma_mat_rel` : file contains in gemma format matrix relatdness used by gemma  (format 1, see manual), matrix must be in same order than fam file. Avoid computation of relatdness by pipeline. 
    * `gemma_multi` : option run gemma by chromosome separately (is not a loco). to increase time computation and used more cpus, 0 no 1 yes [default : 0 (No)]
 *  `boltlmm`: should boltlmm be used? 
    *  see [manual](https://data.broadinstitute.org/alkesgroup/BOLT-LMM/)
    * if SNPs is higher than 950,000, 950,000 SNPs are chosen randomly to build the model (see --modelSnps option in bolt)
    * `covariates_type` / `bolt_covariates_type` : for bolt need to define if covariate is binary (0) or continue (1), a comma-separated list as same order than covariates 
    * `bolt_ld_score_file` : A table of reference LD scores for boltlmm is needed to calibrate the BOLT-LMM statistic (option in boltlmm --LDscoresFile),to choice a specific column in Ld file you can use `bolt_ld_scores_col` option (by default : LDSCORE) if option is not provided --LDscoresUseChip used.
    * `bolt_use_missing_cov` : option to "missing indicator method", by default no activate (0), to activate (1) (--covarUseMissingIndic option in boltlmm), which adds indicator variables demarcating missing status as additional covariates.
    * `bolt_num_cores` if bolt is used set this up to 8
    * `bolt_mem_req` memory required for boltlmm, (default : 6GB)
    * impute2 data in bolt  :
      * bolt_impute2filelist : list of impute2 files, each line contains : `chronumber` `file`, file must be in full pattern
      *`bolt_impute2fidiid` : list of individual in same order than bolt_impute2filelist
 *  `fastlmm`: should fastlmm be used?
    *  see [manual](https://github.com/MicrosoftGenomics/FaST-LMM)
    * `fastlmm_num_cores`: if fastmll is used set this up to 8
    * `fastlmm_mem_req`: memory required for fasttlmm (default : 15GB)
    * `fastlmm_multi` : memory used by fastmll is very big and increase with snp number, option run fastlmm by chromosome, with relatedness matrix computed before with gemma (-gk 1) 
    * `fastlmmc_bin` : should change binary for fastlmmc (default fastlmmc)
and then for all the tests except _gemma_, _boltlmm_ and _fastlmm_, do you want to adjust for multiple testing 
 * plink gwas option :
    * `adjust`: do we want to do explicit testing for Bonferroni correction et al that PLINK odes
    * `mperm`: do you want to test doing permutation testing. If so, how many tests?  By default this is 1000.
 * `fastGWA` :  do an analyse with fastGWA (GCTA), specific option [fastGWA manual](https://cnsgenomics.com/software/gcta/#fastGWA)
   * fastgwa  (if 1 perfom fastGWA, otherwise no [default 0])
   * `fastgwa_type` : [default : --fastGWA-mlm-exact]
   * `fastgwa_memory` : memory for fastgwa and grm [default: 10G] 
   * `fastgwa_cpus` : cpus for fastgw and grm [default: 5]
   * `covariates_type` : similar to `bolt_covariates_type`, give for each covariable type of covariable qualitatif (0) or quantitatif (1), must be len of covariates, if nothing, consider all covariable just as quantitatif covariable  [default ""] 
   *  grm :
    * `gcta_grmfile` : index of file with grm with `gcta_grmfile`.grm.id and `gcta_grmfile`.grm.sp, if extension will not here, grm will build see below
    * to build grm :
     * `grm_nbpart` : nb part to build grm [default : 100]
     * `gcta64_bin` : binary for gcta64 [default : gcta64] 
     * `grm_cutoff` : cutoff value for  grm matrix (using option --make-bK-sparse) [default : 100]
  * `gcta64_bin` : binary for gcta64 [default : gcta64] 

with pipeline, do a GxE interaction with Gemma and Plink, arguments :
  * `gxe` : environmental variables to do gxe analysis with `pheno`, must be coded in 1 and 2 for plink
  * `gemma_gxe` : GxE interation with gemma [default : 0], see  `covariates` to add covariates in gemma models
   * pipeline computed frequencies, N for each group with plink and add to files
  * `plink_gxe` : GxE interation with plink (see option -gxe, in [plink manual](http://zzz.bwh.harvard.edu/plink/anal.shtml#qtgxe)) [default : 0], no covariate could be provided.
   * pipeline computed frequencies, N for each group with plink and add to files, futhermore they add A1 and A2.
   * furthermore BetaGxE and SeGxE computed by pipeline as SeGxE : sqrt(se1^2 + se2^2), BetaGxE : Z\_GXE * SeGxE

For example

```nextflow run assoc/assoc.nf --input_pat raw-GWA-data --chi2 1 --logistic 1 --adjust 1```

analyses the files `raw-GWA-data` bed, bim, fam files and performs a chi2 and logistic regression test, and also does multiple testing correction.

Other flags are:
* `thin`. You can set this to a floating point number in the range (0, 1] and then the PLINK data files are thinned leaving only that proportion of the SNPs. This allows pipeline to be tested with a small proportion of the data This is probably only needed for debugging purposes and usually this should not be be set.
* `chrom`. Only do testing on this chromosome.


# Post-Analysis script 
##1. Computed a p.value par permutation with gemma 
TODO
=======
# 2. Post-Analysis script 

## 2.1. Computed a p.value par permutation with gemma 

## 2.2. Conditional & joint (COJO) analysis of GWAS summary statistic

this section describes a pipeline in devloment, objectives is doing a conditional and joint association using GWAS summary data and gcta
see [cojo](https://cnsgenomics.com/software/gcta/#COJO)

### Installation
need python3, gcta 
tested for singularity image: no
### Running
The pipeline is run: `nextflow run assoc/annot-assoc.nf`

The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * `output_dir` : output directory
  * `output` : output pattern
  * ̀ data` : same option that _assoc/main.nf_, file is optional, used if need select specific individus for gcta,  compute frequencies or N, if mission in `file_gwas`
  * `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
  * `pheno` : optional, header in data, if present select individuals with no missiong individual to keep individuals for computed frequencie or gcta
  * `cut_maf` minor allele frequencies [ default : 0.0001]
  * ̀`file_gwas` : file contains gwas result, if N or frequencies is not available, it is computed with plink file and `data` file, to change format header must be defined :
    * ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
    * `head_freq` : freq header [ default : None], if not present computed with plink, (and data/pheno if present)
    * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
    * `head_rs` : rs header column [default : "SNP"]
    * `head_beta` : beta header colum [default : "BETA"]
    * `head_se`  : column for standard error of beta "SE"
    * `head_A1` : column for A0 :[default : "ALLELE0" ]
    * `head_A2` : column for A1 :[default : "ALLELE1" ]



Cojo parameter :
  * `cojo_wind` :  Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified. [ default : 10000 ]
  * `cojo_actual_geno` : If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. *option to avoid due to a various bug*  [default 0]
  * `cojo_slct` : Perform a stepwise model selection procedure to select independently associated SNPs? 1 : yes 0 : no [default 1]
    * `cojo_p` :  Threshold p-value to declare a genome-wide significant hit. The default value is 5e-8 if not specified. This option is only valid in conjunction with the option `cojo_slct`.
    * `cojo_slct_other` : other option for slct see [manual](https://cnsgenomics.com/software/gcta/#COJO)
  * `cojo_top_snps_chro` :  Perform a stepwise model selection procedure to select a fixed number of independently associated SNPs by chromosome without a p-value threshold.  [integer between 0 and n, to define top snp number. default : 0].
  * `gcta_mem_req`="6GB"


## 2.3. Annotation of position 

This section describes a pipeline in devlopment, objectives is annotation of rs using annotation, locuszoom, and phenotype in function of genotype
>>>>>>> 8933c79753c01b026291ad5777a84455a37dcd3f

### Installation
need locuszoom, _R_ : (ggplot2), python3


## 2.4. Simulation pipeline: `assoc/simul-assoc.nf`

This section describes a pipeline in devlopment, purpose of this pipeline is to estimate false positive and false negative with simulated phenotype, Our script, *assoc/simul-assoc.nf* takes as input PLINK files that have been through quality control and
  * Simulate quantitative phenotypes with [phenosim](https://www.ncbi.nlm.nih.gov/pubmed/21714868) based on genetics data
  * perform a GWAS on  phenotype simulated using gemma, boltlmm.
  * Perform summary statistics.

### Installation
a version of _phenosim_ adapted is already in nextflow binary, write in python2. plink, gemma and bolt must be installed

### Running

The pipeline is run: `nextflow run assoc/simul-assoc.nf`

The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * input, output and script directories: the default is that these are subdirectories of the `work_dir` and there'll seldom be reason to change these;
  * `input_pat` : this typically will be the base name of the PLINK files you want to process (i.e., do not include the file suffix). But you could be put any Unix-style glob here. The workflow will match files in the relevant `input_dir` directory;
  * `num_cores` : cores number used
  * `ph_mem_req` : memory request for phenosim
  *  Simulation option :
     * `phs_nb_sim` : simulation number (default : 5)
     * `phs_quant_trait` :  quantitative trait simulation : 1, qualitative not develop yet (default : 1, -q option in phenosim)
     * Quantitative trait option :
        * `ph_nb_qtl` : number of simulated QTN (default: 2, option -n in phenosim)
        * `ph_list_qtl` : proportion of variance explained by each QTNs, separate the values by commas (default : 0.05 -q in phenosim)
        * `ph_maf_r` :  MAF range for causal markers (upper and lower bound, separated by a comma, no space) (default: 0.05,1.0, -maf_r in phenosim)
        * option to do a linear transformation of phenotype with co factor of external data and normatisation:
           * ph_normalise : perform a normalisation (1) or not 0 (Default)
           * each phenotype i be normalise using newpheno = norm(pheno)+var0i*a+var1i*b+ ... + intercept
           * `ph_cov_norm` : contains coefficients for relation separed by a comma (ex "sex=0.2,age=-0.1)
           * `data` : contains cofactor data for each individuals used to normalise with
           * `ph_cov_range` : normalisation range for initial phenotype
           * `ph_intercept` : intercept
  * Association option :
     * `boltlmm` : 1 perform boltlmm (default 0), see boltlmm option in _assoc/main.nf_
     * `gemma` : 1 perform gemma (default 0)  see gemma option in _assoc/main.nf_
     * `covariates` : covariates to include in model (if ph_normalise is 1)
  * Statistics option :
     * `ph_alpha_lim` : list of alpha used to computed significance (separated by comma)
     * `ph_windows_size` : windows size around position used to simulate phenotype to define if was detected, in bp ex 1000bp in CM ex 0.1CM

### output 
different output is provided :
   * simul folder : contains position used to defined phenotype
   * in boltlmm/gemma folder,  res_boltlmm/gemma.stat  contains summary stat for each alpha:
      * we defined `windows true` as the windows around snp used to build phenotype (size is defined previously)
      * `nsig_simall_alpha` : number significant snp in all windows true
      * `nsig_sim_alpha` :   number windows true where at least one snps is significant
      * `nsig_simaround_alpha` : number significant windows true where one snp is significant and has been excluded snps used to build pheno
      * `nsig_nosim_alpha` : snp significant snp not in windows true
      * `nsnp` : snp number total  in dataset
      * `nsnpsima` : snp number used to build phenotype (see ph_nb_qtl)
   * in boltlmm/gemma/simul/ : contains p.value compute for each simulation

###Note 
  * for phenotype simulation all missing values is discarded and replaced by more frequent allele
  * phenosim use a lot of memory and time, subsample of snp/samples improve times / memory used

