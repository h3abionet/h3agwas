<img src="../../auxfiles/H3ABioNetlogo2.jpg"/>

# H3Agwas Pipeline Version : replication using gwas catalog



## Input :
 * `gwas_file` : summary statistics, need to be described
  * `head_pval`  ["P_BOLT_LMM"]
  * `head_freq`   [""]
  * `head_bp`  "BP"
  * `head_chr` = "CHR"
  * `head_rs` = "SNP"
  * `head_beta`  [""]
  * `head_se` [""]
  * `head_A1` ["ALLELE1"]
  * `head_A0` ["ALLELE0"]
 *  gwas catalog : algorithm download version of gwas catalog from ucsc, separate by comma
  * `pheno` : un phenotype
  * `file_pheno` : file containing list of phenotype extracted from gwas catalog
  * `list_chro` : list chro to analyse, separate by x-y (x to y) and comma [1-22]
 * `wind`
 * `min_win`

## Algorithm :
### list of phenotypes available
one option possible to obtained list of phenotype

"""
nextflow run h3abionet/h3agwas/replication/gwascat/main.nf --justpheno 1 --output_dir gwascatpheno/ --output gwascatpheno -resume -profile slurm
"""

### GWAS catalog :
extract from gwas catalog postiion and chromosome relative to `file_pheno` or `pheno` and `list_chro`. computed heritabilite, beta and se for input
OR transformed using log2 and defined if column contains information higher than 1


### positions replication
algorithm :

* extract from gwas catalog information relative to phenotype defined previously
* extract from gwas file (`gwas_file`) list corresponding to file found in gwas catalog and around each positoins  

### ld replications 
algoritm :
* computed ld between each positions of gwas catalog and summary statistics
* merged 
* keep each positions where p-value < `min pvalue`
* reported each output

#### windows replication in fn
### heritabilities computed





