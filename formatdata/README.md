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

