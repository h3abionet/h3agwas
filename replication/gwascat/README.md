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


### exact replication by positions

algorithm :
* extract from gwas catalog information relative to phenotype defined previously
* extract from gwas file (`gwas_file`) list corresponding to positions found in gwas catalog 
* merge and computed p value adjusted with bonferonni or fdr 
* result in `resultat/exact_rep` folder

### ld replications 
algoritm :
* computed ld between each positions of gwas catalog and position from summary statistics.
* merged positions of summary statistics in ld with gwas catalog  and gwas catalog.
* keep each positions where p-value < `min pvalue` and r2 < `clump_r2`.
* for all positions in LD with one position gwas catalog, min p-value extracted from summary statistics.
some specificity :
* position of gwas catalog cannot be in sumstat file but in genetics file and LD can be computed
* One positions from sumstat can be in ld with more than one position of gwas catalog
* result in `result/ld`

### clump replcation 
algoritm :
* used clump function of  plink with min pvalue and `clump_r2` to defined groups of snps 
* each input of gwas catalog and summary statistics will be merge by groups
* keep p-value clump (min from ld)
* result in `result/clump/`

### windows using kb
algoritm :
* around each position of gwas catalog extract around  `size_win_kb` information from positions keep min p-value
* result in `result/wind`

### windows using ld 

* algoritm windows :
 * run clump function of plink using `clump_r2` and `size_win_kb` using summary statistics
 * keep group of 
 * extracted all positons from summary position between begin, and end used min
 * extracted min p-value by windows
* result in `result/ldwind/wind`

* algoritm wind extended :
 * computed ld between positon  using `clump_r2`, and windows size of `size_win_kb`
 * defined different groups of positions 
 * merged groups where two positions are in two group
 * keep group positons where gwas catalog positon
 * defined windows using chro, min position and max position
 * extracted all positons from summary position between begin, and end used min
 * extracted min p-value by windows 
* result in `result/ldwind/windext`



### heritabilities computed
for summary statistics and gwas catalog computed  using beta, se, af, n
 * (2 x β^2 x MAF x (1−MAF))/((2 x  ^ 2 x MAF x (1−MAF) ) + (se^2 * 2 * N* MAF* (1−MAF))) see [here](://storage.googleapis.com/plos-corpus-prod/10.1371/journal.pone.0120758/1/pone.0120758.s001.pdf?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=wombat-sa%40plos-prod.iam.gserviceaccount.com%2F20210519%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20210519T133314Z&X-Goog-Expires=86400&X-Goog-SignedHeaders=host&X-Goog-Signature=d16ff5bd5266f8caf0656cffd3385dbe10a0341258cd7e2eb82ad325075cc9a7b1aa84d5a9b4b09fe3ff335389c10c7d8c6d1aea893a2ca0f3994dde4c9b002ceab0c79fc0e97b6ad498b998d9a2bcc75cc2be0ceeb3d39baff4b781f1e4f885c08caa844274f7411eef821146519a7e1fe0b0679a67e222217de7d0ce2ab84ea72c6fc08555d49559cb597acfdb2d1c45bb53f2e338e84f6fcc675d2c2d5ca11818711d115b86724e0ec1f2fbd6206705fa3ebd227c0f10d608f25738256c6405f551c81b1603aea753cd9cced1f28e3a249d5de986fece222508826f051900d898ac77f99a8aff97e47ab00572cf9855792eecaccc09c51177771584086d85)
for gwas catalog se defined using IC, as
(Upper - beta) / 1.96

z values defined as beta / se







