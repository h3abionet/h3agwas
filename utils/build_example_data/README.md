<img src="../../../auxfiles/H3ABioNetlogo2.jpg"/>

#  build genotype using 1000 Genotype pipeline: `utils/build_example_data/main.nf`

This workflow has been extensively expanded by Jean-Tristan Brandenburg

## Running

The pipeline is run: `nextflow run utils/build_example_data/main.nassocf`

The key options are:
* `output_dir` : output direction
* `output` : [default "out"]
* `pos_allgeno` : position extracted from all genotype 
* `list_chro`  : list chro extracted fromm 1000 genome

* information relative of gwas catalog to extract:
// * `gwas_cat` : file of gwas catalog, for the moment just [uscs format](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz) is allowed (not inplemented)
 * `gwas_cat_ftp` : file to download gwas catalog  
 * `list_pheno` : list pheno extracted from gwas catalog, split by one comma each pheno (default : "Type 2 diabetes")
* extraction and simulation :
 * `simu_hsq` : see docutmentation of hsq
 * `simu_k` : see docutmentation of gcta 0.01
 * `simu_rep` : repettition
 * `gcta_bin` : gcta for simulation
 * clump position of gwas catalog :
  * `clump_p1` [default :0.0001]
  * `clump_p2` [default :0.01]
  * `clump_r2` [ default : 0.50]
  * `clump_kb` :[ default : 250 ]
* other : 
 * `plk_cpus` : [default 10]
 * `gcta_bin` : [default gcta64]


For example

``` TODO ```


### Installation
need tabix, TODO
tested for singularity image: no
### Running


#  Simulation using phenosim and bed filee: `utils/build_example_data/main.nf`

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

