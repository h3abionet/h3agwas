<img src="auxfiles/H3ABioNetlogo2.jpg"/>


## What's new :
* 2022-07-13 : found example and command line of pipeline in [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)
* 2022-07-13 : add in  [finemapping/finemap_region.nf](finemapping/README.md) , 1000 genome downlowd if there is no genetics data
* 2022-07-13 : add regenie software for association [assoc](assoc/README.md) 
* 2022-07-8 : add plink support for saige [assoc](assoc/README.md) 
* 2022-07-8 : add bgen support for saige,bolt-lmm and gcta [assoc](assoc/README.md) 
* 2022-06-10 : update meta of metasoft, add column for N, Frequencies for each publication in merge result
* 2022-05-04 : preprint manuscript can be found on [biorxiv](https://www.biorxiv.org/content/10.1101/2022.05.02.490206v1)
* 2021-11-04 : association script include saige software
* 2021-10-25 : new script to do a conditional analysis using gemma see [finemapping](finemapping/README.md)
* 2021-09-15 :  add for script in finemapping `finemap_multi.nf`, extract positions using plink, do define independant position in function of minimum p-value and for each independant position apply finemapping [finemapping](finemapping/README.md)
* 2021-03-09 : add plot and merge for estimation of heritabilities
* 2021-03-08 : re organise pipeline :
  * estimation of heritabilities assoc/esth2-assoc.nf move at heritabilites/main.nf
  * do a cojo using gcta : assoc/cojo-assoc.nf move at  finemapping/cojo-assoc.nf
  * meta association assoc/meta-assoc.nf move at meta/meta-assoc.nf
  * meta using mtag association assoc/mtag-assoc.nf move at meta/mtag-assoc.nf
  * permutation assoc/permutation-assoc.nf move at utils/permutation/main.nf
  * simulation of phenotype using phenosim assoc/simul-assoc.nf utils/build_example_data/simul-assoc.nf
* 2021-02-18: add pipeline to build a example data using gwas catalog and 1000 genome [build\_example\_data](utils/build_example_data/README.md)
* 2021-02-16: add report to vcf in plink with analyse of frequencies and score  [formatdata](formatdata/README.md)
* 2021-01-22: create utils folder to add Metasoft binary and utils (server down)
* 2020-12-08: add meta analyse with plink [assoc](assoc/README.md)
* 2020-12-01: add plink GxE, add estimation of beta and se [assoc](assoc/README.md)
* 2020-11-17: add module nf to convert vcf in bgen format [formatdata](formatdata/README.md)
* 2020-07-27: add covariable qualitatif to fastgwa [assoc](assoc/README.md)
* 2020-07-27: News nextflow modules to transform vcf impute format in bimbam[formatdata](formatdata/README.md)
* 2020-06-03: News nextflow modules to transform plink file in vcf file with check allele for imputation[formatdata](formatdata/README.md)
* 2020-05-18: fixed bug in gcta to computed heribilities [assoc](assoc/README.md)
* 2020-03-27: added a modules to convert position between different genome version [formatdata](formatdata/README.md)
* 2020-02-20: support for awsbatch
* 2020-02-20 :  added fastgwa (software gcta) as assoc software  : [assoc](assoc/README.md)
* 2019-10-01 : added in transform data a nextflow script to format output of GWAS with added your own rs, frequencies, N etc...  (usefull for post analysis) : [formatdata](formatdata/README.md)
  * file `formatdata/format_gwasfile.nf`
* 2019/09/19 : added in estimation of heritabilites option for Multiple variance components for boltlmm  [assoc](assoc/README.md)
* 2019/09/17 : added format and analysis by mtag in [assoc](assoc/README.md)
* 2019/09/16 : added two news nextflow files to convert data in [formatdata](formatdata/README.md):
  * `formatdata/vcf_in_plink.nf` : format data in vcf for plink
  * `formatdata/vcf_in_impute2.nf` : extract impute2 data from vcf of sanger
* 2019/09/10 : update estimation of heritability in [assoc](assoc/README.md) to take account for each software when heritabilities can't be computed

