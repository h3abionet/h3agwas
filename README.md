# h3agwas


This is currently in draft

## Background

h3aGWAS is a simple human GWAS analysis workflow originally built at the [Sydney Brenner Institute](http://www.wits.ac.za/academic/research/sbimb/20747/wits_bioinformatics.html) for data quality control (QC) and basic association testing, and later refined and extended by H3ABionet. It uses Nextflow as the basis for workflow managment, and has been dockerised.

* NB: our pipelines are now  under active development. If you would like to do some beta testing of the new pipeline, please checkout the `wits` branch. * If you do so, please ignore the wiki, the documentation is all in the README.


## Documentation 

A
A
A
A
A


Installation, Examples and tutorials for witsGWAS can be found in the wiki

## Features

**QC of Affymetrix array data** (SNP6 raw .CEL files)

  * genotype calling
  * converting birdseed calls to PLINK format

**Sample and SNP QC of PLINK Binaries**

Sample QC tasks checking:

 *  discordant sex information
 *  calculating missingness
 *  heterozygosity scores
 *  relatedness
 *  divergent ancestry 

SNP QC tasks checking:

 * remove duplicates
 * discordant sex information
 * minor allele frequencies
 * SNP missingness
 * differential missingness
 * Hardy Weinberg Equilibrium deviations

**Association testing**

 * Basic PLINK association tests, producing manhattan and qqplots
 * CMH association test - Association analysis, accounting for clusters
 * permutation testing
 * logistic regression
 * emmax association testing



**Running the pipeline**

The pipeline is controlled through the nextflow.config file. All parameters including input files, and parameters.  This can be edited manually 


## Copyright

### Authors

Lerato E. Magosi, Kiran Anmol, Shaun Aron, Rob Clucas, Eugene de Beste, Scott Hazelhurst, Aboyomini Mosaku, Don Armstrong and the Wits Bioinformatics team

### License
witsGWAS is offered under the MIT license. See LICENSE.txt.

### Download

`git clone https://github.com/h3abionet/h3agwas`

### References
