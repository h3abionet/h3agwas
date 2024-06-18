## Prepare data for imputation : `plk_in_vcf_imp.nf`

### argument :
* `input_dir`
* `input_pat`
* `output_dir` : direction of output [default : output]
* file to extract rsinfomation with position :
 * `file_ref_gzip` : must be in gzip example of file used : [here](ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz)
    * `poshead_chro_inforef` psotion of column chromosome in file  [default : 0]
    * `poshead_bp_inforef` : position of column where bp in file [default : 1]
    * `poshead_rs_inforef` : position of column where rs in file  [default : 2]
* `deleted_notref` : deleted position s not found in `file_ref_gzip`
* `reffasta` : fasta reference, if present do control of vcf file :
* checkVCF.py
* bcftools : used plugin of +fixref see `BCFTOOLS_PLUGINS=bcftools/plugins/`
* `michigan_qc` : apply migigan qc [default 0]
 * see : [prepare you data](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/)
 * `dataref_michigan` : data ref used by michigan, if empty dowload [default : ""]
 * `ftp_dataref_michigan` : dl data from michigan [default : ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz]
 * `bin_checkmich` perl script [default : "HRC-1000G-check-bim.pl"]

### General requirement 
* bcftools
* plink
* R
* python
* qctools (v2)
* samtools
* for control of vcf
 * checkVCF.py is present in binary of nextflow pipeline (https://github.com/zhanxw/checkVCF)
 * if michigan qc apply, need data set of michigan and perl script (see [here](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/))

### Example 

see [h3agwas-example github](https://github.com/h3abionet/h3agwas-examples)
```
 nextflow run h3abionet/h3agwas/formatdata/vcf_in_plink.nf --file_listvcf utils/listvcf --output_pat  kgp_imputed --output_dir plink_imputed/   --reffasta utils_data/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz  -profile singularity
```

# format vcf
## input 
* format vcf in various output :
 * plink
 * bgen
 * bimbam
* if any input vcf from previous 
 * `vcf_folder_zip` : zip folder contained all folder as same pipeline if password used `vcf_folder_zip_password`
 * `vcf` : one vcf file 
 * `vcf_list`:file contained list of vcf
