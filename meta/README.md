<img src="../auxfiles/H3ABioNetlogo2.jpg"/> 
#  Meta analyse pipeline: `assoc/main.nf`

## Mtag analysis

### reference 

### parameters
TODO

  * `file_gwas` : one ore more one file gwas of differents phenotype
    * ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
    * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
    * `head_rs` : rs header column [default : "SNP"]
    * `head_beta` : beta header colum [default : "BETA"]
    * `head_se`  : column for standard error of beta "SE"
    * `head_A1` : column for A0 :[default : "ALLELE0" ]
    * `head_A2` : column for A0 :[default : "ALLELE2" ]
    * `head_freq` : freq header [ default : A1Freq],
    * `head_n`: N header, used just for ldsc, if not present, `Nind` must be initialize.
  * if n not initialise :
    * used plink file to computed each position with n :
      * `input_pat` : input pattern of plink file
      * `input_dir` : input dir of plink file
    * list_n : need to be implemented

