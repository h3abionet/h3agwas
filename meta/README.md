<img src="../helperfiles/H3ABioNetlogo2.jpg"/> 
#  Meta analyse pipeline: `assoc/main.nf`

# 3. MetaAnalysis pipeline : `assoc/meta-assoc.nf`

This section describes a pipeline in devlopment, purpose of this pipeline is to do a meta analysis with a various format files.Our script, *meta-assoc.nf* takes as input various GWAS results files and `rsid` to do a metanalysis with METAL, GWAMA and Metasoft

### Installation
need python3, METAL (last version : https://github.com/statgen/METAL), GWAMA, MR-MEGA and MetaSoft (one version is available on utils/bin/)

### Option

The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * `input`, `output` and script directories: the default is that these are subdirectories of the `work_dir` and there'll seldom be reason to change these;
  * `output_dir` = "all"
  * meta analysis option :
     * `metal` : 1 perform metal (default 0)
     * `gwama` : 1 perform gwama (default 0)
     * `plink` : 1 perform perform meta analyse in plink(default 0)
     * `metasoft` : 1 perform metasoft(default 0)
       * `metasoft_pvalue_table` : for metasoft need files :  _HanEskinPvalueTable.txt_
     * `mrmega` : 1 perform MR-MEGA (default 0)
  * `file_config`
     * describe all informations for each gwas result used for meta analysis
     * file is comma separated (csv), each line is to describe one file
     * header of config file is : rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,IsRefFile
       * `rsID` : column name for rsID in gwas file
       * `Chro` : column name for Chro in gwas file
       * `Pos` : column name for Pos in gwas file
       * `A1` :  column name for reference allele in gwas file
       * `A2` :  column name for alternative allele in gwas file
       * `Beta` :  column name for B values in gwas file
       * `Se` :  column name for sterr values in gwas file
       * `N` : column name for size in gwas file
       * `freqA1` : column name for freqA1 or maf in gwas file
       * `direction` : column name of strand for association -/+  in gwas file
       * `Imputed` :  column name of imputed or not for position in gwas file
       * `NCount` : column name to add a column N at your file with value in the column
       * `Sep` : what separator is in gwas file :
         * you could use characters as ; . : but to avoid some trouble you can use :
           * COM : for comma
           * TAB : for tabulation
           * WHI : for white space
       * `File` : gwas file with full path
       * if one of the column is missing in your GWAS file, replace by _NA_
  * optional option :
     * memorie usage :
       * plink_mem_req : [20GB]
       * gwama_mem_req : gwama memories [20GB]
       * metasoft_mem_req : metasoft memories ["20G"]
       * ma_mem_req : request for extraction of data, change format and plot of manhathan ["10G"]
       * mrmega_mem_req : mr mega memorie ["20GB"]
     * cpu memorie :
        * max_plink_cores : [default 4]
        * other used 1 cpus
     * binaries :
       * `metal_bin` : binarie for metal (default : _metal_ )
       * `gwama_bin` :  binarie for gwama ( default : _GWAMA__ )
       * `metasoft_bin` : binarie for java of metasoft ( default _Metasoft.jar_)
       * `mrmega_bin` : binarie for java of metasoft ( default _Metasoft.jar_)
       * `plink_bin` : binarie for java of metasoft ( default _Metasoft.jar_)
     * options software :
       * `ma_metasoft_opt` : append other option in metasoft command line(default : null)
       * `ma_genomic_cont` : use a genomic_control use in METAL and GWAMA(default, 0)
       * `ma_inv_var_weigth`: do a invert variance weight usefull for metal (default, 0)
       * `ma_overlap_sample`: do you have sample overlap? used by metal(default, 0)
       * `ma_random_effect` : do mixed model (default 1)
       * `ma_mrmega_pc` : how many pcs used for mrmega (default : 4)
       * `ma_mrmega_opt` : append other option in MR-MEGA command line (default : null)
    * `us_rs` : if you want chromosome and position are replaced using rs (warning you need to be sure that one chromosome position has same rs in each file), [default 0, yes : 1], otherwise they will used chromosome and position to replaced rs
### specificity 
#### MR-MEGA
MR-MEGA need chromosomes, positions and N (sample number) for each position, so in pipeline referent file (in file_config, 1 in IsRefFile) must be have chromosome and poosition

### option comparaison software

| | options | 
| Software | `ma_genomic_cont` | `ma_inv_var_weigth` | `ma_overlap_sample` | `ma_random_effect` |
| --- | --- | --- | --- | --- |
| descriptionus| genomic control | invert variance weight | sample overlap | Random Effect |
| default | 0 | 0 | 0 | 0 |
| metal | yes | yes  | yes | no |
| gwama | yes | default  no | yes | 
| Mr Mega| yes | no | no | no |  
| plink | no | yes | no | no |
| Metasoft | no | no | no | default |
| --- | --- | --- | --- | --- |

1 'weighted-z' requests weighted Z-score-based p-values (as computed by the Abecasis Lab's METAL software)

#### Ressource

| Software | Manuals | References |
| --- | --- | --- |
| Metal | [here](https://genome.sph.umich.edu/wiki/METAL_Documentation) |  [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922887/)|
| Gwama | [here](https://genomics.ut.ee/en/tools) | [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-288)| 
| Mr Mega| [here](https://genomics.ut.ee/en/tools)|  [here](https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/1471-2105-11-288.pdf)|
| plink | [here](https://zzz.bwh.harvard.edu/plink/metaanal.shtml) | [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/)|
| Metasof | [here](genetics.cs.ucla.edu/meta) | [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3146723/) |
| --- | --- | --- |

### References
[Comparison of Two Meta-Analysis Methods: Inverse-Variance-Weighted Average and Weighted Sum of Z-Scores](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5287121/)

### Example 
* data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)
* a csv file need to described each input, contains header for each file

```
echo "rsID,Chro,Pos,A1,A2,Beta,Se,Pval,N,freqA1,direction,Imputed,Sep,File,Ncount" > utils/input_meta.csv
for File in `ls data/summarystat/*.gemma|grep -v ".all"`
do
echo "rs,chr,ps,allele0,allele1,beta,se,p_wald,NA,af,NA,NA,TAB,$File,2500" >>  utils/input_meta.csv
done
```

* input :
  * user can choose software that he want to run : metal (`--metal 1`), gwama (`--gwama 1`), metasoft (` --metasoft 1`) MrMega (`--mrmega 1`) and plink (`--plink 1`)

```
nextflow run ~/Travail/git/h3agwas/meta/meta-assoc.nf   --metal 1 --gwama 1 --metasoft 1 --mrmega 1 --plink  1  --file_config utils/input_meta.csv -resume -profile slurmSingularity --output_dir meta
```



## Mtag analysis

### reference 

### parameters

* `file_gwas` : one ore more one file gwas of differents phenotype
  *  ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
  * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
  * `head_rs` : rs header column [default : "SNP"]
  * `head_beta` : beta header colum [default : "BETA"]
  * `head_se`  : column for standard error of beta "SE"
  * `head_A1` : column for A0 :[default : "ALLELE0" ]
  * `head_A2` : column for A0 :[default : "ALLELE2" ]
  * `head_freq` : freq header [ default : A1Freq],
  * `head_n`: N header, used just for ldsc, if not present, `Nind` must be initialize.
* if n not initialise :
  * used plink file to computed each position with n :
  * `input_pat` : input pattern of plink file
  * `input_dir` : input dir of plink file
  * list_n : need to be implemented


## Example
* data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)
* multi-trait analysis of GWAS (MTAG), a method for joint analysis of summary statistics from genome-wide association studies (GWAS) of different traits, possibly from overlapping samples.
* input :
 * list of summary statistic `file_gwas` and header from gwas file: `-head_[name]`
 * also you can give nformation relative to : ` --input_dir data/imputed/ --input_pat imput_data --pheno pheno_qt1,pheno_qt2 --data data/pheno/pheno_test.all `, can add N value to each summary statistic

```
nextflow h3abionet/h3agwas/meta/mtag-assoc.nf --head_freq af --head_pval p_wald --head_bp ps --head_chr chr --head_rs rs --head_beta beta --head_se se --head_A1 allele1 --head_A2 allele0 --input_dir data/imputed/ --input_pat imput_data --file_gwas  data/summarystat/all_pheno.gemma,data/summarystat/all_phenoq2.gemma --pheno pheno_qt1,pheno_qt2 --data data/pheno/pheno_test.all -resume   -profile slurmSingularity
```

