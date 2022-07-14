<img src="../../auxfiles/H3ABioNetlogo2.jpg"/>
#  Annotation pipeline : `annotation/main.nf`

This section describes a pipeline in devlopment, purpose of this pipeline is to extract, plot of specific position from gwas result, phenotype and genotype, and annotate positoin

## Main option
The key options are:
* `work_dir` : the directory in which you will run the workflow. This will typically be the h3agwas directory which you cloned
* `output_dir` : output directory
* `output` : output pattern
* `data` : same option that _assoc/main.nf_, your file of phenotype and covariable
 *  `pheno` : phenotype, must be in data file
 * `cov` : covariable, if not null must be in data
* `input_pat` : the base of set of PLINK bed,bim and fam files (this should only match one);
* `input_dir` : the directory of set of PLINK bed,bim and fam files (this should only match one);
* `file_gwas` : file contains gwas result, if N or frequencies is not available, it is computed with plink file and data file, to change format header must be defined :
 * `head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
 * `head_freq` : freq header [ default : None], if not present computed with plink, (and data/pheno if present)
 * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
 * `head_rs` : rs header column [default : "SNP"]
 * `head_beta` : beta header colum [default : "BETA"]
 * `head_se` : column for standard error of beta "SE"
 * `head_A1` : column for A0 [default : "ALLELE0" ]
 * `head_A2` : column for A1 [default : "ALLELE1" ]
* `list_rs` :
 * list of rs where want an analyse, same name of rs must be present in gwas file
* Annotation :
 * `file_annotation` :
  *  file contains list file by chromosome with information for annotation
  * todescribe
 * `info_annotation` :
  * give information for each annotation 
  * to described
* locuszoom :
 * `loczm_buid` : locus zoom build (default hg19)
 * `loczm_pop`  : population of locus zoom
 * `loczm_source` : source of locus zoom
 * `loczm_gwascat` :  [defaut : none ]
 * `loczm_bin` : binary of locus zoom, used to defined also database 
## Installation 
 * locuszoom, R : (ggplot2), python3

## Example 
* Data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)

### warning : 
 aws : locus zoom doesn't work, pipeline on AWS doesn't work, you must install your own locus zoom

### command lines
input :
 * plink file
 * summary statistisc
 * phenotype and data file
 * rs from gwas (one or more, separate with comma)
 * Annotation used annovar file, if you don't `list_file_annot` and `info_file_annot`, data will be downloaded

```
nextflow run  ~/Travail/git/h3agwas/utils/annotation/main.nf --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --input_dir data/imputed/  --input_pat imput_data --file_gwas data/summarystat/all_pheno.gemma  --output_dir annotation --list_rs "2:45832137:A:G,1:117539108:G:T" --data data/pheno/pheno_test.all --pheno pheno_qt1  -resume  -profile slurmSingularity --loczm_bin  "/dataE/AWIGenGWAS/shared/ResultGWAS/Ressource/locuszoom/bin/locuszoom"
```


