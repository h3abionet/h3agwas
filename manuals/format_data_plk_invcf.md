# Prepare data for imputation : `plk_in_vcf_imp.nf`

## Arguments :
### input 
* `input_dir`
* `input_pat`
* `bfile`
* file to extract rs with position :
 * `file_ref_gzip` : must be in gzip example of file used : [here](ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz)
    * `poshead_chro_inforef` psotion of column chromosome in file  [default : 0]
    * `poshead_bp_inforef` : position of column where bp in file [default : 1]
    * `poshead_rs_inforef` : position of column where rs in file  [default : 2]

* `deleted_notref` : deleted position s not found in `file_ref_gzip`

### output 
* `output_dir` : direction of output [default : output]
* `output` : direction of output [default : output]

### software
* bcftools : used plugin of +fixref see `BCFTOOLS_PLUGINS=bcftools/plugins/`
* `fasta` : fasta reference, if present do control of vcf file 


* checkVCF.py
 * see : [prepare you data](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/)

### General requirement 
* bcftools
* plink
* R
* python
* samtools
* for control of vcf
 * checkVCF.py is present in binary of nextflow pipeline (https://github.com/zhanxw/checkVCF)
 * if michigan qc apply, need data set of michigan and perl script (see [here](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/))

### Example 

see [h3agwas-example github](https://github.com/h3abionet/h3agwas-examples)
```
 nextflow run h3abionet/h3agwas/formatdata/vcf_in_plink.nf --file_listvcf utils/listvcf --output_pat  kgp_imputed --output_dir plink_imputed/   --fasta utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz  -profile singularity
```

