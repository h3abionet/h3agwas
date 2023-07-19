<img src="../../../auxfiles/H3ABioNetlogo2.jpg"/>

#  build genotype using 1000 Genotype pipeline and gcta: `utils/build_example_data/main.nf`

## Running

The pipeline is run: `nextflow run utils/build_example_data/main.nf`

The key options are:
* `output_dir` : output direction
* `output` : [default "out"]
* `pos_allgeno` : file contains chromosome and position that will be extracted from vcf file, each lines is chromosome and position [need, no default]
* `list_chro`  : list chro extracted fromm 1000 genome [1-22,X : from  1 to 22 and X chromosome]
* ftp1000genome : download 1000 genome for build dataset no : 0 yes :1 [default :1]
* `list_vcf` : file contains list of vcf if ftp1000genome is no,

* information relative of gwas catalog to extract:
 * `gwas_cat` : file of gwas catalog, for the moment just [uscs format](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz) is allowed (not inplemented)
 * `gwas_cat_ftp` : file to download gwas catalog  
 * `list_pheno` : list pheno extracted from gwas catalog, split by one comma each pheno (default : "Type 2 diabetes")
* extraction and simulation :
 * `simu_hsq` : variance explained of disease by genetics [default 0.3]
 * `simu_k` :   Specify the disease prevalence for binary phenotype [default 0.1]
 * `simu_rep` : repetition number using same snps and heriability effect [default 10]
 * `used_effect` :  two behavious for effect, used gwas catalog effect (option 1) or simulated value using a normal law  (option 0) [default : 1]
 * clump position of gwas catalog :
  * `clump_p1` [default :0.0001]
  * `clump_p2` [default :0.01]
  * `clump_r2` [ default : 0.50]
  * `clump_kb` :[ default : 250 ]
* other : 
 * `plk_cpus` : cpus number used for plink and gcta[default 10]
 * `gcta_bin` : binary for gcta [default gcta64]


### Example
#### created input file 
```
awk '{print $1"\t"$4}' data/array_plk/array.bim  > utils/list_posarray
```
#### command line 

```
nextflow run ~/Travail/git/h3agwas//utils/build_example_data/main.nf -profile slurmSingularity   --pos_allgeno utils/list_posarray -resume --nb_snp 3 --output_dir simul_gcta_main
```


### Installation
need tabix, TODO
tested for singularity image: yes


#  Simulation using gcta and plink file : `utils/build_example_data/simul-assoc_gcta.nf`
Pipeline used same argument than previous pipeline without download 1000 genomes 


## Example 
* Data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)

```
nextflow run ~/Travail/git/h3agwas/utils/build_example_data/simul-assoc_gcta.nf -profile slurmSingularity  --input_dir data/imputed/  --input_pat  imput_data --output_dir simul_gcta
```




#  Simulation using phenosim and plink file : `utils/build_example_data/simul-assoc_phenosim.nf`

This section describes a pipeline in devlopment, purpose of this pipeline is to estimate false positive and false negative with simulated phenotype, Our script, *assoc/simul-assoc.nf* takes as input PLINK files that have been through quality control and
  * Simulate quantitative phenotypes with [phenosim](https://www.ncbi.nlm.nih.gov/pubmed/21714868) based on genetics data
  * perform a GWAS on  phenotype simulated using gemma, boltlmm.
  * Perform summary statistics.

### Installation
a version of _phenosim_ adapted is already in nextflow binary, write in python2. plink, gemma and bolt must be installed

### Option

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

### Output 
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
### Example :
* Data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)
```
nextflow run ~/Travail/git/h3agwas/utils/build_example_data/simul-assoc_phenosim.nf -profile slurmSingularity  --ph_normalise 0 --input_dir data/imputed/ --input_pat  imput_data --gemma 1
```



