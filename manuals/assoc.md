
# Main option
## Phenotype options
* `data` : data files [default previous workflow]
* `phenotype` : phenotypes
* `phenotype_type` : phenotypes type binary or continuous
* `covariates_type` :  type of covariates
* `pheno_tr_fct` : transform phenotypes using a functions [default : '']
* `add_pcs` : number pcs to add as covariates [default 0]
* `pheno_residuals` : your phenotypes in residuals using covariates [default : 0]
* `phenores_tr_fct` transform residuals using function [default : '']
 * function to transform your variabes (`pheno_tr_fct`, `phenores_tr_fct`):
   * any R function to transform data
   * log, log10... 
   * invnorm  : inverse normal rank transformation
## Genetics

## option association :
* `assoc` : performed association 
 * `saige`
 * `regenie`
* other option
 * `gxe` : variable to performed a GXE [default : '']

## genetics  
* `loco` [default 1]
* `plink_indep_pairwise` [100 20 0.1]
* `cut_maf_rel` [0.01]

## memory and cpu association

## specific option association
### Saige
* `saige_bin_fitmodel` ["step1_fitNULLGLMM.R"]
* `saige_bin_spatest` ["step2_SPAtests.R"]
* `saige_otheropt_step1` ['']
* `saige_otheropt_step2` ['']
* `saige_otheropt`  ['']
* `saige_impute_method` ["best_guess"]

## 

| | option effect | gxe |
| gxe  | --gxe |  yes | 
| loco | --loco |  yes  |

## G 

