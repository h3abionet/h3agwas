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
* `genetic_map_file` : genetic map used in boltlmm 
* build relatdness and GRM, you have different way to buil matrix of relatdness for boltlmm, gemma, fast gwa or fastlmm
  * by default all snps are used and for boltlmm 9,500,000 are shuffled 
  * you can give a rs contains id or rs `file_rs_buildrelat` 
  * you can used `sample_snps_rel` (default 0), will used plink to sample snps and 
  *  `file_rs_buildrelat` : file with rs list (one by lines) to build genetics models (relatdness), for gemma `-snps` for boltlmm `--modelSnps`


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


