# h3agwas


This is currently in draft

## Background

h3aGWAS is a simple human GWAS analysis workflow originally built at the [Sydney Brenner Institute](http://www.wits.ac.za/academic/research/sbimb/20747/wits_bioinformatics.html) for data quality control (QC) and basic association testing, and later refined and extended by H3ABionet. It uses Nextflow as the basis for workflow managment, and has been dockerised.

The original version of the h3aGWAS was published in June 2017 with minor updates and bug fixes through the rest of the year. Based on experience with large data sets, the pipelines were considerably revised with additional features, reporting and a slightly different workflow.  

There is one feature of the original workflow that has been omitted. Version 1 supported parallel GWAS anaysis of different data files in one Nextflow run. This has been removed. Although, not unuseful, this feature complicated the implementation and made expansion more difficulkt and also this capacity can be simulated easily at the operating system level (see the Wiki).

We have moved all scripts from Python 2 to Python 3, so you will need to have Python 3 installed. We are working on moving away from the use of R and Perl to simplify installation.

*This is a development version -- it has been tested on real data but the containerised workflow has not yet been impemtented*

## Documentation 

Installation, Examples and tutorials for witsGWAS can be found in the wiki

## Features

**QC of Illumina array data**

  * Converting from an Illumina topbottom format report into PLINK

**Sample and SNP QC of PLINK Binaries**

Sample QC tasks checking:

 *  discordant sex information
 *  calculating missingness
 *  heterozygosity scores
 *  relatedness
 *  divergent ancestry 

SNP QC tasks checking:
 * batch reports
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
