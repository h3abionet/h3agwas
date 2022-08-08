<img src="../helperfiles/H3ABioNetlogo2.jpg"/>


# Pipelines Format data, plink in vcf, vcf in various format, summary statistics and change of build

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
  * give a report with analyse of frequencie and score
### arguments :
* `file_listvcf` : file contains each bgzip vcf files to merge, one by line [default : none]
* `min_scoreinfo` : what score info minimum do you want : [default :0.6]
* `output_pat` : pattern of output for bed final file [default : out]
* `output_dir` : directory of output : [default : plink] 
* `score_imp` : header of score imputation [default : INFO], for score info depend of software used on software of imptuation
 * PBWT : INFO
 * 
* `do_stat` : by default true make stat using frequencies and score 
 * statfreq_vcf : 
  * pattern used in Info to computed frequencies  ([default : "%AN %AC" with AN total and AC alternative number]) 
  * can be two value NAll NAlt, where frequencies computed as Nalt/NAll
  * can be one value frequencies
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

###  Example 
* data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)


## Script : `vcf_in_bgen.nf` / `vcf_in_bgen_merge.nf` `vcf_in_bgen_merge_chro.nf`
### three srategy of conversion :
 * `vcf_in_bgen.nf` : convert each vcf that was filter in bgen, without merge with one process 
 * `vcf_in_bgen_merge.nf` : 
   * vcf are filter and merge and merge and format in bgen
 * `vcf_in_bgen_merge_chro.nf` :
  * vcf are filter
  * vcf are format in bgen 
  * bgen are merge
### Requirement :
* plink, bcftools, bash, qctools, nextflow
* singularity / dockers image : no test yet
### what script done :
* Intial data : format of Sanger imputation format vcf file
* select for each chromosome on quality of imputation : min info
* convert each vcf in impute2 format used by boltlmmm
* output for each chromosome is basename of initial files with .impute2.gz

### arguments :
* `file_listvcf` : file contains each bgzip vcf files to merge, one by line [default : none]
* `min_scoreinfo` : what score info minimum do you want : [default :0.6]
* `output_dir` : directory of output : [default : bgen] 
* `qctoolsv2_bin` : bin file for qctools 
* `genotype_field` : genotype field to transform [degault : GP]
* `bgen_type` : bgen type see [qctool manual] (https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/alphabetical_options.html) :
 * default bgen (other  : "bgen_v1.2", "bgen_v1.1")
* `other_option`  : other option to give to qctools

for instance for bolt lmm format bgen must be :
```
~/nextflow ~/Travail/git/h3agwas/formatdata/vcf_in_bgen_merge.nf --file_listvcf listvcf --output_pat  exampledata2_imp --output_dir ./ -profile slurmSingularity -resume --bgen_type bgen_v1.2 
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

## Pipeline : convert position hg38 in hg19
```
nextflow run convert_posversiongenome.nf
```
### what is doing?
* if no file give download gwas catalog
* extract positions of interest 
* used rs to search position see args `file_ref_gzip`
* used crossmap to defined position s not found previously and strand : see `bin_crossmap` and `data_crossmap`
* return file with new position
### arguments
* `output_dir` : direction of output [default : output]
* `output `: output : [default : out]
* `file_toconvert` : file to convert if empty download gwas catalog
  * `link_gwas_cat`  : link to download gwas catalog [default : https://www.ebi.ac.uk/gwas/api/search/downloads/alternative ]
  * `head_rs` : head rs to file to convert [default SNPS (gwas catalog)] 
  * `head_bp` : head bp to file to convert [default SNPS (gwas catalog)] 
  * `head_chro` : head bp to file to convert [default SNPS (gwas catalog)] 
  * `sep` : separator used TAB, SPACE, "," [default TAB] (not allowed : ;)

* file to extract rsinfomation with position :
 * `file_ref_gzip` : must be in gzip example of file used : [here](ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/All_20180423.vcf.gz)
    * `poshead_chro_inforef` psotion of column chromosome in file  [default : 0]
    * `poshead_bp_inforef` : position of column where bp in file [default : 1]
    * `poshead_rs_inforef` : position of column where rs in file  [default : 2]

* `bin_crossmap` : crossmap [default ~/.local/bin/CrossMap.py]
* `data_crossmap` : data to convert [default : "" ]
  * if no argument will be download : 
   * hg38 in hg19 : `link_data_crossmap` (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz)

###output ;
* `{out}.tsv` : final file 
* `{out}.multi.tsv` : more that one position have been found 
* `{out}.detail.tsv` :  file before cleaning
* `{out}.notfound.tsv` : fileswhere position not found
* folder `datai` : contains files contains files download
* folder `datatmp` : contains temporary file (extract of rs file)
### installation :
R : library
pip3.6 install CrossMap --user
pip3.6 install numpy==1.16.1 --user
chmod +x ~/.local/bin/CrossMap.py

## Pipeline : `vcf_in_bimbam.nf`
transform vcf in bimbam format after filters for quality.
###arguments
* `file_listvcf` : file contains each bgzip vcf files to merge, one by line [default : none]
* `min_scoreinfo` : what score info minimum do you want : [default :0.6]
* `output_dir` : directory of output : [default : impute2]
* `genotype_field` : genotype field in vcf file [default : GP]
* `qctoolsv2_bin`  : qctools v2 binary [default :qctool_v2]
* `bcftools_bin` : bcftools bin [default : bcftools]


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
 

# General requirement 
* bcftools
* plink 
* R
* python
* qctools (v2)
* samtools
* for control of vcf 
 * checkVCF.py is present in binary of nextflow pipeline (https://github.com/zhanxw/checkVCF)

#Example 

see [h3agwas-example github](https://github.com/h3abionet/h3agwas-examples)
