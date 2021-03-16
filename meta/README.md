<img src="../auxfiles/H3ABioNetlogo2.jpg"/> 
#  Meta analyse pipeline: `assoc/main.nf`

# 3. MetaAnalysis pipeline : `assoc/meta-assoc.nf`

This section describes a pipeline in devlopment, purpose of this pipeline is to do a meta analysis with a various format files.Our script, *meta-assoc.nf* takes as input various GWAS results files and `rsid` to do a metanalysis with METAL, GWAMA and Metasoft

### Installation
need python3, METAL, GWAMA, MR-MEGA and MetaSoft

### Running
The pipeline is run: `nextflow run meta-assoc.nf`

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
       * `Sep` : what separator is in gwas file :
         * you could use characters as ; . : but to avoid some trouble you can use :
           * COM : for comma
           * TAB : for tabulation
           * WHI : for white space
       * `File` : gwas file with full path
       * `IsRefFile` : you need to define a reference file to define what rs should be considered in other files
       * if one of the column is missing in your GWAS file, replace by _NA_
  * optional option :
     * binaries :
     * binaries :
       * `metal_bin` : binarie for metal (default : _metal_ )
       * `gwama_bin` :  binarie for gwam ( default : _GWAMA__ )
       * `metasoft_bin` : binarie for java of metasoft ( default _Metasoft.jar_)
       * `mrmega_bin` : binarie for java of metasoft ( default _Metasoft.jar_)
       * `plink_bin` : binarie for java of metasoft ( default _Metasoft.jar_)
     * options software :
       * `ma_metasoft_opt` : append other option in metasoft command line(default : null)
       * `ma_genomic_cont` : use a genomic_control use in METAL and GWAMA(default, 0)
       * `ma_inv_var_weigth`: do a invert variance weight usefull for metal (default, 0)
       * `ma_random_effect` : do mixed model (default 1)
       * `ma_mrmega_pc` : how many pcs used for mrmega (default : 4)
       * `ma_mrmega_opt` : append other option in MR-MEGA command line (default : null)
### specificity 
#### MR-MEGA
MR-MEGA need chromosomes, positions and N (sample number) for each position, so in pipeline referent file (in file_config, 1 in IsRefFile) must be have chromosome and poosition



## Mtag analysis

### reference 

### parameters
TODO

  * `file_gwas` : one ore more one file gwas of differents phenotype
    * ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
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

