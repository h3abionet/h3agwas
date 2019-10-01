<img src="../aux/H3ABioNetlogo2.jpg"/>

# Pipeline : format data

## Script : `vcf_in_plink.nf`
### Requirement :
  * python3, plink, bcftools, bash, nextflow
  * singularity / dockers image : no test yet
### what script done :
  * select for each chromosome on quality of imputation : min info
  * convert each vcf in plink
  * rename duplicate rs or . by chro
  * added cm in bim files if file `genetic_map` give in argumen, 
  * merge all chromosome in plink format
### arguments :
* `file_listvcf` : file contains each bgzip vcf files to merge, one by line [default : none]
* `min_scoreinfo` : what score info minimum do you want : [default :0.6]
* `output_pat` : pattern of output for bed final file [default : out]
* `output_dir` : directory of output : [default : plink] 
* `genetic_maps` : genetics maps to added map in bim file, if not provided, map doesn't added in bim, must be not compressed :
 * file for [hg19](https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg19.txt.gz) 
 * file for [hg17](https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg17.txt.gz)
 * file for [hg18](https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg18.txt.gz)
 * file for [hg38](https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz)

```
chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)
1 55550 0 0
1 568322 0 0
1 568527 0 0
```


### example 

```
ls  dirvcf/chr*.vcf.gz > listfileawigen
nextflow run h3abionet/h3agwas/formatdata/vcf_in_plink.nf --file_listvcf listfileawigen -resume -profile slurm --output_pat awigen  --genetic_maps $FileCM --plink_mem_req 6GB -r hackathon
```

## Script : `vcf_in_impute2.nf`
### Requirement :
* plink, bcftools, bash, nextflow
* singularity / dockers image : no test yet
### what script done :
* Intiial data : format of Sanger imputation format vcf file
* select for each chromosome on quality of imputation : min info
* convert each vcf in impute2 format used by boltlmmm
* output for each chromosome is basename of initial files with .impute2.gz

### arguments :
* `file_listvcf` : file contains each bgzip vcf files to merge, one by line [default : none]
* `min_scoreinfo` : what score info minimum do you want : [default :0.6]
* `output_dir` : directory of output : [default : impute2] 

## format a gwas file `format_gwasfile.nf`
### Requirement :
* plink, bash, nextflow, python3 (library : panda)
* singularity / dockers image : no test yet

### what script done :
* initial data of gwas format transform in other files
* search rs on file to added new rs at each position (if not found add chro:pos)
* added N and frequencies values if need and plink file gave
* Change header, separator... etc

### arguments :
* `file_gwas` : one file gwas :
  * intial header of your file :
    * `head_pval` [optional]
    * `head_freq` [optional]
    * `head_bp`  
    * `head_chr` 
    * `head_beta` [optional]
    * `head_se` [optional]
    * `head_A1` [optional]
    * `head_A2` [optional]
    * `head_N`  [optional]
  * `sep` separator default space or tab, [optional], for comma : COM, tabulation : TAB and space "WHI"
* header of your output :
  * if not initialise : using output of your initial files
  * `headnew_pval` [optional]
  * `headnew_freq` [optional]
  * `headnew_bp`  [optional]
  * `headnew_chr` [optional]
  * `headnew_beta` [optional]
  * `headnew_se` [optional]
  * `headnew_A1` [optional]
  * `headnew_A2` [optional]
  * `headnew_N`  [optional]
* file to extract rsinfomation with position :
  * `file_ref_gzip` : must be in gzip example of file used : [here](ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz)
    * `poshead_chro_inforef` psotion of column chromosome in file  [default : 0]
    * `poshead_bp_inforef` : position of column where bp in file [default : 1]
    * `poshead_rs_inforef` : position of column where rs in file  [default : 2]

* others option :
  * added some N and frequencie in gwas file :
    * used plink information to compute freq and N and added in gwas file if `head_N` and/or `head_freq` not intialise
    * `input_dir` : plink directory
    * `input_pat` : plink basename
  * `mem_req` : memory request for processes>

