<img src="../aux/H3ABioNetlogo2.jpg"/>

#  finemapping pipeline: `finemapping/main.nf`

This workflow has been extensively expanded by Jean-Tristan Brandenburg

The purpose of this pipeline is to perform a initial analysis of finemapping 

Our script, *finemapping* takes as input PLINK files, gwas file

## Running

The pipeline is run: `nextflow run finemapping`

The key options are:
* `input_dir`, `output_dir`: where input and output goes to and comes from;
* `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
* `file_gwas` : file contains summary stat of gwas :
 * ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
 * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
 * `head_rs` : rs header column [default : "SNP"]
 * `head_beta` : beta header colum [default : "BETA"]
 * `head_se`  : column for standard error of beta "SE"
 * `head_A1` : column for A0 :[default : "ALLELE0" ]
 * `head_A2` : column for A0 :[default : "ALLELE2" ]
 * `head_freq` : freq header [ default : A1Freq],
 * `head_chr` : freq header [ default : A1Freq],
 * `head_bp` : freq header [ default : A1Freq],

* `n_pop` : number individuals for data set pop [default : 10000]

* `gwas_cat` : file of gwas catalog for plot 
  * `headgc_chr` : chromosome header of gwas catalog
  * `headgc_bp` : position head of gwas catalog

* **binary** :
 * `finemap_bin` : finame binary 
 * `caviar_bin` : finame binary 
 * `caviar_bin` : finame binary 
 * `n_causal_snp` : for finemapping cond and xxx number causal max snps used
* Cojo parameter :
  * `gcta_mem_req`="6GB"
  * `cojo_wind` :  Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified. [ default : 10000 ] (need to be implemented)
  * `cojo_actual_geno` : If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. *option to avoid due to a various bug*  [default 0] (need to be implemented)
  * `cojo_slct_other` : other option for slct see [manual](https://cnsgenomics.com/software/gcta/#COJO) )(need to be implemented)


### Installation
need locuszoom, _R_ : (ggplot2), python3

For example

TODO

analyses the files `raw-GWA-data` bed, bim, fam files and performs a chi2 and logistic regression test, and also does multiple testing correction.



