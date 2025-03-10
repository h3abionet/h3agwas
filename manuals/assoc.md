
# Main option
## phenotype options
* `data`
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
## option test :
* `assoc` : performed association 
 * `saige`
 * `regenie`
* other option
 * `gxe` : variable to performed a GXE [default : '']
## genetics  
 * `loco`

## spe
### Saige
 * `saige_otheropt_step1`

## 

| | option effect | gxe |
| gxe  | --gxe |  yes | 
| loco | --loco |  yes  |

## G 

