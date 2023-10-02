<img src="../helperfiles/H3ABioNetlogo2.jpg"/>

# Finemaping 

5 scripts to perform fine-mapping :
  * run fine-mapping at one specific region `finemapping/finemap_region.nf` 
  * run fine-mapping on full summary statistics `finemapping/main.nf` 
  * Stepwise model selection procedure to select independently associated SNPs with  [cojo](https://cnsgenomics.com/software/gcta/#COJO) from gcta `finemapping/cojo.nf` 
  * Conditional analysis using gemma `finemapping/cond-assoc.nf`
  * Conditional analysis using gcta and summary statistics `finemapping/cond-gcta.nf`
 
##  Finemapping a specific region using pipeline: `finemapping/finemap_region.nf`

The purpose of this pipeline is to perform a initial analysis of finemapping 

Our script, *finemapping* takes as input PLINK files, gwas file

### some limits 
* pipeline discarded positions duplicated from genetics file and summary statistics
### Running

The pipeline is run: `nextflow run finemapping/finemap_region.nf`

### options

The key options are:
* finemapping :
  * `n_causal_snp` : for finemapping causal snp number
  * `chro` : chromosome to apply chromosome
  * `begin_seq` : begin sequence to apply sequence
  * `end_seq` : end sequence to apply sequence
  * `prob_cred_set` :  prob of credible set [default : 0.95], for FINEMAP software
  * `used_pval_z` : build beta and se using z or pvale [defaul 0:no]
* `threshold_p` : treshold to used [default : 5E10-8]
* genetics data :
  * 2 way are possible, data already in plink format (if possible same data that association has been done) or 1000 genome will be download
  * genetics data user :
    * `input_dir`, `output_dir`: where input and output goes to and comes from;
    * `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
* 1000 genome :
  * `ftp_vcf`
  * `other_cpus_req`
* `file_gwas` : file contains summary stat of gwas :
  * `head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
  * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
  * `head_rs` : rs header column [default : "SNP"]
  * `head_beta` : beta header colum [default : "BETA"]
  * `head_se`  : column for standard error of beta "SE"
  * `head_A1` : column for allele 1 :[default : "ALLELE1" ]
  * `head_A2` : column for allele 0 or 2 :[default : "ALLELE0" ]
  * `head_freq` : freq header [ default : ""],
  * `head_chr` : freq header [ default :  ""],
  * `head_bp` : freq header [ default : ""],
* `n_pop` : number individuals for data set pop [default : 0]
* `gwas_cat` : file of gwas catalog for plot 
  * `headgc_chr` : chromosome header of gwas catalog
  * `headgc_bp` : position head of gwas catalog
  * `list_phenogc` :  from gwas catalog what pheno used (separated by a comma)
  * `file_phenogc` :  from gwas catalog what pheno used inside a file (separated by a comma)

* cojo parameter for Stepwise model selection procedure to select independently associated SNPs:
  * `gcta_mem_req`="6GB"
  * `cojo_wind` :  Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified. [ default : 10000 ] (need to be implemented)
  * `cojo_actual_geno` : If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. *option to avoid due to a various bug*  [default 0]
  * `cojo_slct_other` : other option for slct see [manual](https://cnsgenomics.com/software/gcta/#COJO) )(need to be implemented)
* annotations parameter :
 * `paintor_fileannot)`  : file contains annotation (see paintor manual)

* **binary** :
  * `finemap_bin` : software binary 
  * `paintor_bin` : software binary 
  * `plink_bin` : software binary 
  * `gcta_bin` : software binary 
* mem and cpu :
  * `max_plink_cores` LD / plink [default 4]
  * `plink_mem_req`   LD / plink     [default : "5GB"]
  * `gcta_mem_req` [default :"6GB"]
  * `caviar_mem_req` [default :"40GB"]
  * `paintor_mem_req` [default :"20GB"]
  * `fm_mem_req ` [default : "20G"]
  * `plink_mem_req` [default :"6GB"]


### Installation
need locuszoom, _R_ : (ggplot2), python3, finemap, paintor, gcta, plink, gwama

### Example 
* Data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)
```
nextflow run  h3abioneth3agwas/finemapping/finemap_region.nf  --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --list_phenogc "Type 2 diabetes" --input_dir  data/imputed/  --input_pat imput_data --file_gwas data/summarystat/all_pheno.gemma  --output_dir finemapping_pheno1_wind --output finemapping_pheno1 -resume  -profile slurmSingularity --begin_seq 112178657 --end_seq 113178657 --chro 10
```

##  finemapping pipeline automated with selection of lead snps using plink : `finemapping/main.nf`

### algorithm :
 * using clump plink to defined independant locus
 * for each snps apply algorithms from pipeline `finemapping/finemap_region.nf`

### options

The key options are:
* finemapping :
  * `n_causal_snp` : for finemapping causal snp number
  * `prob_cred_set` :  prob of credible set [default : 0.95], for FINEMAP software
  * `used_pval_z` : build pvalue using z (need to check) [defaul 0:no]
  * `threshold_p` : treshold to used [default : 5E10-8]
  * `threshold_p2` :  Secondary significance threshold for clumped SNPs see [clump](https://zzz.bwh.harvard.edu/plink/clump.shtml)
  * `size_wind_kb` : size windows used for clump and around each independant snps
  * `clump_r2` : ld used for [clump](https://zzz.bwh.harvard.edu/plink/clump.shtml)
  * `cut_maf` minor allele frequencies [ default : 0.01]
* `data` :
  * `data` : optional used individual from data to clean plink file`
  * `pheno` : phenotype associate to data
  * `covariates` : covariate must be in data
* plink :
  * `input_dir`, `output_dir`: where input and output goes to and comes from;
  * `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
* `file_gwas` : file contains summary stat of gwas :
  * `head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
  * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
  * `head_rs` : rs header column [default : "SNP"]
  * `head_beta` : beta header colum [default : "BETA"]
  * `head_se`  : column for standard error of beta "SE"
  * `head_A1` : column for allele 1 :[default : "ALLELE1" ]
  * `head_A2` : column for allele 0 or 2 :[default : "ALLELE0" ]
  * `head_freq` : freq header [ default : TODO],
  * `head_chr` : freq header [ default : TODO],
  * `head_bp` : freq header [ default : TODO],
* `n_pop` : number individuals for data set pop [default : 0]
* `gwas_cat` : file of gwas catalog for plot
  * `headgc_chr` : chromosome header of gwas catalog
  * `headgc_bp` : position head of gwas catalog
  * `list_phenogc` :  from gwas catalog what pheno used (separated by a comma)
  * `file_phenogc` :  from gwas catalog what pheno used inside a file (separated by a comma)
* cojo parameter (Stepwise model selection procedure to select independently associated SNPs.) :
  * `gcta_mem_req`="6GB"
  * `cojo_wind` :  Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified. [ default : 10000 ]
  * `cojo_actual_geno` : If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. *option to avoid due to a various bug*  [default 0] 
  * `cojo_slct_other` : other option for slct see [manual](https://cnsgenomics.com/software/gcta/#COJO) 
* mem and cpu :
  * `max_plink_cores` LD / plink [default 4]
  * `plink_mem_req`   LD / plink     [default : "5GB"]
  * `gcta_mem_req` [default :"6GB"]
  * `caviar_mem_req` [default :"40GB"]
  * `paintor_mem_req` [default :"20GB"]
  * `fm_mem_req ` [default : "20G"]
  * `plink_mem_req` [default :"6GB"]
* annotations parameter :
  * `paintor_fileannot`  : file contains annotation (see paintor manual)

### Installation
need locuszoom, _R_ : ggplot2, python3, finemap, paintor, gcta, plink

### Output 
* $output/clump/clump_output.clumped result of clump by plink 
* $output/gwascat : file of gwas catalog format
* $output/data : gene dowmloaded
* $output/fm/$chr\_$bp : each positions where finemapping had been apply :
  * caviarbf :result of caviarbf
  * cojo_gcta :  result a stepwise model selection procedure to select independently associated SNPs. showing the LD correlations between the SNPs. using gcta
  * fm\_cond : finemap software using option using `--cond` 
  * fm\_sss  : finemap software using option `--sss`
  * paintor : paintor output 
  * $output/$chr\_$bp/Finemapping\_.all.out : result merge of all result
    * information relative at position : rsid    chromosome      position        allele1 allele2 rsid    chromosome      position        allele1 allele2
    * result of association : beta, se, p (if used `used_pval_z = 1`, beta and se had been computed using p-value, otherwise orignel)
    * gcta result, givin your independant : bj   bJ_se    pJ    LD_r    logpJ 
    * value of your association rsid	chromosome	position	allele1	allele2	maf	beta	se
    * `IsSig` : position significant
    * `is_cred` : is in credible interval set using finemap with `--sss` 
    * `cred_$chr_$bp_cond.cred` : is in credible of finemap for finemap software with option `--cond`
    * `cred_$chr_$bp_sss.cred` : is in credible of finemap for finemap software with option `--sss`
    * `finemap` sss result : `prob_fm_sss`, `log10bf_fm_sss`,	`mean_fm_sss`, `sd_fm_sss`,`mean_incl_fm_sss`,	`sd_incl_fm_sss`



### Example

* Data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)

```
nextflow run  h3abioneth3agwas/finemapping/main.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --list_phenogc "Type 2 diabetes" --input_dir  data/imputed/  --input_pat imput_data --file_gwas data/summarystat/all_pheno.gemma  --output_dir finemapping_pheno1 --output finemapping_pheno1 -resume  -profile slurmSingularity
```




## Conditional & joint (COJO) analysis of GWAS summary statistic
this section describes a pipeline in devlopment, objectives is doing a Stepwise model selection procedure to select independently associated SNPs on full summary statistics
see [cojo](https://cnsgenomics.com/software/gcta/#COJO)

### Installation
need python3, gcta
tested for singularity image: yes

### Running
The pipeline is run: `nextflow run finemapping/annot-assoc.nf`

The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * `output_dir` : output directory
  * `output` : output pattern
  * `data` : same option that _assoc/main.nf_, file is optional, used if need select specific individus for gcta,  compute frequencies or N, if mission in `file_gwas`
  * plink :
   * `input_dir`, `output_dir`: where input and output goes to and comes from;
   * `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
  * `pheno` : optional, header in data, if present select individuals with no missiong individual to keep individuals for computed frequencie or gcta
  * `cut_maf` minor allele frequencies [ default : 0.0001]
  * ̀`file_gwas` : file contains gwas result, if N or frequencies is not available, it is computed with plink file and `data` file, to change format header must be defined :
    * `head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
    * `head_freq` : freq header [ default : None], if not present computed with plink, (and data/pheno if present)
    * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
    * `head_rs` : rs header column [default : "SNP"]
    * `head_beta` : beta header colum [default : "BETA"]
    * `head_se`  : column for standard error of beta "SE"
    * `head_A1` : column for A0 :[default : "ALLELE0" ]
    * `head_A2` : column for A1 :[default : "ALLELE1" ]

* Cojo parameter :
  * `cojo_wind` :  Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified. [ default : 10000 ]
  * `cojo_actual_geno` : If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. *option to avoid due to a various bug*  [default 0]
  * `cojo_slct` : Perform a stepwise model selection procedure to select independently associated SNPs? 1 : yes 0 : no [default 1]
    * `cojo_p` :  Threshold p-value to declare a genome-wide significant hit. The default value is 5e-8 if not specified. This option is only valid in conjunction with the option `cojo_slct`.
    * `cojo_slct_other` : other option for slct see [manual](https://cnsgenomics.com/software/gcta/#COJO)
  * `cojo_top_snps_chro` :  Perform a stepwise model selection procedure to select a fixed number of independently associated SNPs by chromosome without a p-value threshold.  [integer between 0 and n, to define top snp number. default : 0].
  * `gcta_mem_req`="6GB"

### Example 
* Data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)

```
nextflow run   h3abionet/h3agwas/finemapping/cojo-assoc.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --input_dir data/imputed/  --input_pat imput_data  --output_dir cojo --data data/pheno/pheno_test.all --pheno pheno_qt1 --file_gwas data/summarystat/all_pheno.gemma  -resume   -profile slurmSingularity
```

## Conditional analysis of GWAS using gemma

`cond-assoc.nf`

### Installation
need python3, gemma, R
tested for singularity image: yes


### algorithm
* pipeline used :
  * plink file and phenotype file 
  * `pos_cond` : list position to be conditionned, transform in 0 (homozygote A1), 0.5 (heterozygote) and 1 (homozygote A2) and used as covariable
  * `chro_cond` : chromosome to be conditionned 
  * `pos_ref` : position de reference where will be extracted positions of interrest 
* pipeline will run :
  * for each `pos_cond`, run gemma on the region using as covariable `pos_cond`
  * `pos_cond`, run gemma on the region using as covariable all `pos_cond` (call merge)
  * run gemma on the region using just phenotype
* output :
  * each gemma output 
  * plot of LD between ref and cond
  * plot of p-value comparison of initial without conditional and conditional 


### Example 
* data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)
* to perform a conditional analysis, using gemma where positions is used as covariable of the phenotype and check pos_ref (`--pos_ref`) is link or indepependant  , argument same than gwas for gemma where you need to add :
  *  pipeline will run a raw gwas using phenotype and covariable and after performed gwas using genotypes of each position from `pos_cond` as covariable
  * `chro_cond` : chro where `pos_cond` and `pos_ref`
  * `pos_ref` : pos will be tested for independance or not to `pos_cond`
  * `pos_cond` : list of position as conditional to verify indpendance with `pos_ref`, include as covariable
* pipeline will compute ld between your positions ref and cond and produce figure


## Conditional analysis of GWAS using gcta

need python3, gcta, R
tested for singularity image: yes

### algorithm
pipeline will extract positions `pos_cond`, `pos_ref` from summary statistics and genetics data and apply options `--cojo-cond` of [gcta](https://yanglab.westlake.edu.cn/software/gcta/#COJO) . 
### main options
* pipeline used :
 * plink file, summary statistics and optional file for individual to clean
 * `pos_cond` : list position to be conditionned
 * `chro_cond` : chromosome to be conditionned
 * `pos_ref` : 1 position de reference where will be extracted positions of interrest
 * `cojo_actual_geno`
 * plink :
  * `input_dir`, `output_dir`: where input and output goes to and comes from;
  * `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
 * `data` : same option that _assoc/main.nf_, file is optional, used if need select specific individus for gcta,  compute frequencies or N, if mission in `file_gwas`
  * `phenotype`
 * `pheno` : optional, header in data, if present select individuals with no missiong individual to keep individuals for computed frequencie or gcta
 * `cut_maf` minor allele frequencies [ default : 0.0001]
 * ̀`file_gwas` : file contains gwas result, if N or frequencies is not available, it is computed with plink file and `data` file, to change format header must be defined :
  * `head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
  * `head_freq` : freq header [ default : None], if not present computed with plink, (and data/pheno if present)
  * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
  * `head_rs` : rs header column [default : "SNP"]
  * `head_beta` : beta header colum [default : "BETA"]
  * `head_se`  : column for standard error of beta "SE"
  * `head_A1` : column for A0 :[default : "ALLELE0" ]
  * `head_A2` : column for A1 :[default : "ALLELE1" ]


