<img src="../../auxfiles/H3ABioNetlogo2.jpg"/>

#  Association pipeline: `annot/main.nf`

## 4. Simulation pipeline: `assoc/simul-assoc.nf`

This section describes a pipeline in devlopment, purpose of this pipeline is to estimate false positive and false negative with simulated phenotype, Our script, *assoc/simul-assoc.nf* takes as input PLINK files that have been through quality control and
  * Simulate quantitative phenotypes with [phenosim](https://www.ncbi.nlm.nih.gov/pubmed/21714868) based on genetics data
  * perform a GWAS on  phenotype simulated using gemma, boltlmm.
  * Perform summary statistics.

### Installation
a version of _phenosim_ adapted is already in nextflow binary, write in python2. plink, gemma and bolt must be installed

