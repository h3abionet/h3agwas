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
  * `list_pheno` : list of phenotype extracted from gwas catalog, if null (all value extracted)
  * `list_chro` : list chro to analyse, separate by x-y (x to y) and comma [1-22]

## Algorithm :
### GWAS catalog :
extract from gwas catalog postiion and chromosome relative to `list_pheno` and `list_chro`. computed heritabilite, beta and se for input
OR transformed using log2 and defined if column contains information higher than 1

### position analyse 
extract from gwas catalog information relative to pheno

### heritabilities computed




