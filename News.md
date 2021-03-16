<img src="auxfiles/H3ABioNetlogo2.jpg"/>


## What's new :
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

