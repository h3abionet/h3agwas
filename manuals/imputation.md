# Imputation modules

## Input
* `vcf` : from previous process
* `build_genome` : ["37"]

* `eagle_genetic_map` :  
 * `ftp_eagle_genetic_map` : ""
 * `ftp_eagle_genetic_map` = [""]
 * `ftp_eagle_genetic_map_38` = ["https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"]
 * `ftp_eagle_genetic_map_19` = ["https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz"]
* `fasta`

* `imp_listchro` :   list of Chro to impute  [default ALL]

* `reference panel` : 
  * `params.imp_ref_name`
  * `params.ref_m3vcf`
    * see [ressource](https://share.sph.umich.edu/minimac3/)
  * ` params.imp_ref_vcf`
* pre - qc :
 * `imp_min_ac`  [default 2]
 * `imp_min_alleles` [default : 2]
 * `imp_max_alleles` [default :2]
* split :
  * `imp_chunk_size`
  * `imp_chunk` : specific chunk separated by a comma ([default ""])
  * `imp_buffer_size` : ([default "100000"]
  * `imp_minRatio` see --minRatio of minimac4 '0.01'
*  
 * `cut_maf_postimp`
 * `impute_info_cutoff`

