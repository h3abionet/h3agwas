<img src="../../auxfiles/H3ABioNetlogo2.jpg"/>

# H3Agwas Pipeline Version : replication using gwas catalog
script using gwas catalog to compared gwas summary statistics to compared gwas summary statistics, using positon, ld, clump or windows based


## Input :
 * `gwas_file` : summary statistics, need to be described
  * `head_pval` : header for pvalue of gwas file ["P_BOLT_LMM"]
  * `head_freq`  :  header for frequencies of gwas file [""]
  * `head_bp` " header of pos of gwas file ["BP"
  * `head_chr` : header of of gwas file "CHR"
  * `head_rs` : header of rsid of gwas file "SNP"
  * `head_beta`  [""]
  * `head_se` [""]
  * `head_A1` ["ALLELE1"]
  * `head_A0` ["ALLELE0"]
 *  gwas catalog : algorithm download version of gwas catalog from ucsc, separate by comma
  * `pheno` : un phenotype
  * `file_pheno` : file containing list of phenotype extracted from gwas catalog
  * `list_chro` : list chro to analyse, separate by x-y (x to y) and comma [1-22]
 * `clump_r2` : min pvalue used for clump and ld 
 * `size_win_kb` : size in kb to analyse used in clump, ld and windows analyse
 * `min_pval_clump` :  p value for clump and defined significant ld algoritm[0.001]
 * `threshold_pval_gwascat` : threshold of pval for gwas catalog (not implemented)

## Algorithm :
### list of phenotypes available
to obtain list of phenotype without run algorinm
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
* keep each positions where p-value < `min pvalue` and r2 < `min r2`
* reported each output merged by gwas catalog input
* position of gwas catalog cannot be in dataset and will be computed distance with other dataset
### clump replcation 
algoritm :
* used clump function of  plink using min pvalue and min r2 with windows to defined clump
* each input of gwas catalog and summary statistics will be merge by windows clump

### windows
algoritm :
* for each position of gwas catalog withh extract around  `size_win_kb` information and position of gwas

### ld clump replication :
algoritm :
 * computed ld between positon in 

#### windows replication in fn
### heritabilities computed





