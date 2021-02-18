<img src="../../../auxfiles/H3ABioNetlogo2.jpg"/>

#  build genotype using 1000 Genotype pipeline: `utils/build_example_data/main.nf`

This workflow has been extensively expanded by Jean-Tristan Brandenburg

## Running

The pipeline is run: `nextflow run utils/build_example_data/main.nassocf`

The key options are:
* `output_dir` : output direction
* `output` : [default "out"]
* `pos_allgeno` : position extracted from all genotype 
* `list_chro`  : list chro extracted fromm 1000 genome

* information relative of gwas catalog to extract:
// * `gwas_cat` : file of gwas catalog, for the moment just [uscs format](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz) is allowed (not inplemented)
 * `gwas_cat_ftp` : file to download gwas catalog  
 * `list_pheno` : list pheno extracted from gwas catalog, split by one comma each pheno (default : "Type 2 diabetes")
* extraction and simulation :
 * `simu_hsq` : see docutmentation of hsq
 * `simu_k` : see docutmentation of gcta 0.01
 * `simu_rep` : repettition
 * `gcta_bin` : gcta for simulation
 * clump position of gwas catalog :
  * `clump_p1` [default :0.0001]
  * `clump_p2` [default :0.01]
  * `clump_r2` [ default : 0.50]
  * `clump_kb` :[ default : 250 ]
* other : 
 * `plk_cpus` : [default 10]
 * `gcta_bin` : [default gcta64]


For example

``` TODO ```





### Installation
need tabix, TODO
tested for singularity image: no
### Running

