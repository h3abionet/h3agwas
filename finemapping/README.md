<img src="../aux/H3ABioNetlogo2.jpg"/>

#  finemapping pipeline: `finemapping/main.nf`

This workflow has been extensively expanded by Jean-Tristan Brandenburg

The purpose of this pipeline is to perform a initial analysis of finemapping 

Our script, *finemapping* takes as input PLINK files, gwas file
## to do
* instalation of programmes  caviarbf, finemapping....
* figures for annotation need included at figures
* caviarbf : need to integrate annotation
* problem with output of option  `--cond` of finemapping
## Running

The pipeline is run: `nextflow run finemapping`

## options

The key options are:
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
* `n_pop` : number individuals for data set pop [default : 10000]
* `gwas_cat` : file of gwas catalog for plot 
  * `headgc_chr` : chromosome header of gwas catalog
  * `headgc_bp` : position head of gwas catalog
* cojo parameter :
  * `gcta_mem_req`="6GB"
  * `cojo_wind` :  Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified. [ default : 10000 ] (need to be implemented)
  * `cojo_actual_geno` : If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. *option to avoid due to a various bug*  [default 0] (need to be implemented)
  * `cojo_slct_other` : other option for slct see [manual](https://cnsgenomics.com/software/gcta/#COJO) )(need to be implemented)
* annotations parameter :
 * Todo
* **binary** :
 * `finemap_bin` : finame binary 
 * `paintor_bin` : finame binary 
 * `plink_bin` : finame binary 
 * `gcta_bin` : finame binary 
 * `n_causal_snp` : for finemapping cond and xxx number causal max snps used

### Installation
need locuszoom, _R_ : (ggplot2), python3

##For example

TODO

## Conditional & joint (COJO) analysis of GWAS summary statistic

this section describes a pipeline in devloment, objectives is doing a conditional and joint association using GWAS summary data and gcta
see [cojo](https://cnsgenomics.com/software/gcta/#COJO)

### Installation
need python3, gcta
tested for singularity image: no

### Running
The pipeline is run: `nextflow run assoc/annot-assoc.nf`

The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * `output_dir` : output directory
  * `output` : output pattern
  * ̀ data` : same option that _assoc/main.nf_, file is optional, used if need select specific individus for gcta,  compute frequencies or N, if mission in `file_gwas`
  * `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
  * `pheno` : optional, header in data, if present select individuals with no missiong individual to keep individuals for computed frequencie or gcta
  * `cut_maf` minor allele frequencies [ default : 0.0001]
  * ̀`file_gwas` : file contains gwas result, if N or frequencies is not available, it is computed with plink file and `data` file, to change format header must be defined :
    * ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
    * `head_freq` : freq header [ default : None], if not present computed with plink, (and data/pheno if present)
    * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
    * `head_rs` : rs header column [default : "SNP"]
    * `head_beta` : beta header colum [default : "BETA"]
    * `head_se`  : column for standard error of beta "SE"
    * `head_A1` : column for A0 :[default : "ALLELE0" ]
    * `head_A2` : column for A1 :[default : "ALLELE1" ]

Cojo parameter :
  * `cojo_wind` :  Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified. [ default : 10000 ]
  * `cojo_actual_geno` : If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. *option to avoid due to a various bug*  [default 0]
  * `cojo_slct` : Perform a stepwise model selection procedure to select independently associated SNPs? 1 : yes 0 : no [default 1]
    * `cojo_p` :  Threshold p-value to declare a genome-wide significant hit. The default value is 5e-8 if not specified. This option is only valid in conjunction with the option `cojo_slct`.
    * `cojo_slct_other` : other option for slct see [manual](https://cnsgenomics.com/software/gcta/#COJO)
  * `cojo_top_snps_chro` :  Perform a stepwise model selection procedure to select a fixed number of independently associated SNPs by chromosome without a p-value threshold.  [integer between 0 and n, to define top snp number. default : 0].
  * `gcta_mem_req`="6GB"




