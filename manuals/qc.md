<img src="../helperfiles/H3ABioNetlogo2.jpg"/>
# Quality control  : 3 mains options
 * deleted and check tecnical duplicate `qc_dup`
 * Quality Control `qc`
 * Clean frequency using an external dataset from  michigan `qc_michigan`

# Clean and check technical duplicate `qc_dup = 1` 

This section describe a way to manage technical duplicate, algoritms is to check if duplicate are correct and deleted putative errors.

## Algorithms
Algorithms :
 * genetics file is clean using data individual and column FID,IID (see `col_fidid`)
 * autosomal relatdness computed between all technical duplicated (see `col_fidiid` and `col_newfidid`)
 * autosomal relatdness computed between all individual where `pi_hat` > `pi_hat_dup`
 * computed quality of duplicate
 * deleted duplicates where `pi_hat` between duplicate  less than `pi_hat_dup` (duplicated are not real duplicated)
 * deleted duplicate and non individual where value `pi_hat` > `pi_hat_dup` (duplicate had been exchange during )
 * for good duplicates where `pi_hat` > `pi_hat_dup`, select best quality of duplicate, deleted others.
 * after cleaning, rename by `col_newfidid` genetic and data (other column in the file will be keep)

## General option
* `qc_dup` : check if duplicate [default 0]
* input : 
 * `data` : phenotype file
  * `col_fidid` : column contained FID and IID same than bfile, 2 columns separated by comma [default : FID,IID]
  * `col_newfidid` : 
    * ID to defined duplicate 
    * same ID in two `individuals` are technical duplicate. 
    * After cleaning of duplicate, initial FID, IID will be replace by `col_newfidiid` in phenotype and genotype
    * one or two heading separated by comma [default : IID_update]
 * genetic :
  * `bfile`
  * `input_dir` and `input_pat`
 * other option :
  * `pi_hat_dup` [default : 0.7]

## Output :
 * report (in work)
 * bfile and data can be used by `qc`, `qc_michigan`, etc..

## example 


## need 
* plink, R
 
# Quality Controls : `qc = 1`

This section describes the various ways in which the pipeline can be run and various options. Usually options are specified in the _nextflow.config_ file (or which ever file you use). However, you can also pass parameters to the Nextflow script on the command-line. Parameters on the command line over-ride any parameters specified in the config file. If you are a first time user of the pipeline you should read the discussion of the `nextflow.config` file in the parent README.


The main pipeline is the PLINK QC pipeline. It takes as input PLINK bed,bim,fam files and performs quality control on  the data according to the parameters specified in the config file.

The Nextflow script file is called *qc.nf*. This could be
called, for example, by running `nextflow run qc`.

The output of the QC is a set of PLINK files that can be used for GWAS, as well as PDF report that describes the QC steps.

### Restrictions
This version has been run on real data sets and works. However, not all cases have been thoroughly tested. In particular
* it is not robust when X chromosome data is not available
* the reporting assumes you want to do batch/site analysis. If you don't the code works but the report may look a bit odd with some figures repeated.

## 1. Overview of the workflow

The QC process consists of:

* removing duplicate markers;
* indentifying indviduals for whom there is discordant sex information;
* removing individuals with too high missingness or excessive heterozygosity;
* detecting whether there are any related individuals and removing enough to ensure that there are not related pairs;
* removing SNPs with too low MAF, or too high missingness, or anomalous HWE, or SNPs where there is a high differential missingness between cases and controls;
* a PCA of the resultant data is computed;
* a detailed report of the QC process is done.

## 2. Input/Output :

### Input 


Users will run the pipeline giving as input PLINK 1.9 bed, bim and fam files.  The key Nextflow parameters to set are:

* Plink and phenotype file Input from a previous part pipeline `qc_dup` 


or

* plink input :
 * `bfile` : basename of plink file, or see `input_pat` and `input_dir`
 * `input_pat` : this typically will be the base name of the PLINK files you want to process (i.e., do not include the file suffix). But you could be put any Unix-style glob here. The workflow will match files in the relevant `input_dir` directory;

*   `batch`: if you want to do QC at a batch level, you need to specify a file with the batch information. This should be a standard PLINK phenotype file (with labels as the first line). If you specify "false" or 0, then no batch-analysis is done. Typically batch information relates to the batches in which samples were genotyped is not intrinsic to the data (e.g. you genotype the first 2500 samples that are available).
* `data` (previously `phenotype`) If you are doing batch analysis you may wish to show how different sub-groups perform in QC with respect to the batch. You will then specify a PLINK-style phenotype file (with labels as the first name).  For example, if you have a multi-site project, you may choose to use the site information as a phenotype. Other possibilities are sex and self-identified group. If you specify "false" or 0, no categorisation will be done. [default : '']
 * `pheno_qc` (previously `pheno_col`) is the column label of the column in  the phenotype file which should be used.
 *  `case_control`  

### Output 
* `output_dir`: the directory which the output should go to. The default is `output`.
* `output`: the base name of the output files. 

## 3. Additional QC Parameters

The following parameters control QC

*  `sexinfo_available`: `true` or `false`. If we don't have sex information then we cannot do the check for discordant genotype. Note that it does not make sense (and is an error) to have sexinfo_available set to true when there is no X-chromosme data in the file;

*  `f_lo_male` and `f_hi_female`. Discordant sex genotype is done on the X-chromosome using the non-recombining parts. F, the in-breeding coefficient of the X-chromosome is computed. If F is above `f_lo_male`, the individual is predicted to be male; if F is below `f_hi_female`, the individual is predicted to be female. Anyone in between is flagged. These cut-off values are arbitrary and especially in large samples you are likely to find a range of F values. However, a large number of discrepant values probably indicates a sample mishandle error.  The PLINK default values (0.8 and 0.2) are the default parameters of the pipeline.
*  `cut_het_high`: What is the maximum allowable heterozygosity for individualsl;
*  `cut_het_low`: minimum
*  `cut_maf `: the minimum minor allele frequency a SNP must have to be included
*  `cut_diff_miss `: allowable differential missingness between cases and controls;
*  `cut_geno`: maximum allowable per-SNP mssingness
*   `cut_mind`: maximum allowable per-individual missingness
*   `cut_hwe`: minimum allowable per-SNP Hardy-Weinberg Equilibrium p-value 
*   `pi_hat`:  maximum allowable relatedness
*   `remove_on_bp`: the first step in the pipeline is to remove duplicate SNPs. There are two ways of detecting duplicates. First, if SNPs have duplicate names (column 1 -- numbering from 0 -- of the bim file). We always remove duplicate SNPs based on this since PLINK gets very upset otherwise. Second, if they are at the same chromosome and base position. If this variable is set to 1, then duplicates based on chromosome or base position are removed too.
 *   `batch_col`: the column label of the file to be used.

 * default `--data` or `--phenotype` if `case_control_col` is intialise  
 * This is the name of a PLINK-style phenotype file with labels in the first line. This is a compulsory parameter. The QC process uses the case/control status of individuals. A principal component analysis is done. We do not expect typically overall that there will be difference between cases and controls. The PC-analysis tests that this is so. Of course, you need to apply your mind to the result as YMMV. If your study has several case/control categories, choose an appropriate one that will give insight. If you only have continuous measures (e.g., BMI), then discretise and make an artificial case-control category. Remember, this is for QC purposes not to find interesting biology. 
  * `case_control_col`: this is the label of the column.
* `high_ld_regions_fname`: this is optional -- it is a list of regions which are in very high LD -- and are exclude when checking for relationships (https://www.cog-genomics.org/plink/1.9/filter#mrange_id).  Provide either absolute file path or relative to where you are running. In a previous version this was relative to input_dir, which is not right.
See [https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)](https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)) for a discussion.

### Chromosome options  (new in version 4)
* `autosome_plink` : how to defined autosome of plink as option [default : --chr 1-22,25 ]
* `chrxx_plink` : number of chr xx in plink format([23])
* `chrxy_plink` : number of chr xy in plink format([25])
* X chromosome filters :
 * `build_genome` : build genome used to split X and XY [default : hg19]
 * `cut_maf_xfemale` [0.01]
 * `cut_maf_xmale` [0.01]
 * `cut_miss_xfemale` [0.02] 
 * `cut_miss_xmale` [0.02] 
 * `cut_diffmiss_x` [0.02]
* Y chromosome filters :
 * `cut_maf_y` : [0.01]
 * `cut_miss_y` : [0.02]

## filtering on frequency :
* we used michigan pipeline appraoche to fitler on frequency, parameters are 
* `input_dir` :
* `input_pat` : 
* `bfile`
* `michigan_qc` : apply migigan qc [default 0]
 * see : [prepare you data](https://imputationserver.readthedocs.io/en/latest/prepare-your-data/)
 * `dataref_michigan` : data ref used by michigan, if empty dowload [default : ""]
 * `ftp_dataref_michigan` : dl data from michigan [default : ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz]
 * `bin_checkmich` perl script [default : "HRC-1000G-check-bim.pl"]

Several of the above parameters make reference to a phenotype file. Of course, these can be to the same phenotype file, but probably using different columns.

### Filtering by GenCall 10 score

If your sample sheet contains a column with the GC10 score (it must be called `10%_GC_Score`), you can filter out all individuals who have a GC10 score below a specified value. 

Please read (the Illumina explanation)[https://www.illumina.com/Documents/products/technotes/technote_gencall_data_analysis_software.pdf] if you are not clear about this.

To do this you need to give a file with the GC10 score (usually it will be in the sample sheet from the sequencing centre). This can either be an Excel or CSV (comma-separated) file.  The header line must contain a column `10%_GC_Score` and a column called `Institute Sample Label` which matches the IDs found in your PLINK file  (an alternative is that there are columns `Sample Label` and `Well` which when concatenated with an underscore give you the ID). The following parameters are relevant

* `samplesheet` : Give the name of the file here. Put 0 (the default) if no samplesheet and then this filtering is not done.
* `gc10`: this is the gc10 score that will be used as the cut-off. The default is 0.4, which is may be too low, but you should definitely think about it.
* `idpat`: Naming conventions differ from sequencing centre to sequencing centre. You may be lucky and the ID in the "Institute Sample Label` matches exactly the FID, IID columns in yoour PLINK data, but more often than not, there is some mangling. For example, it's often the case that the "Institute Sample Label" value contains the sample plate, well _and_ your study ID (e.g. WP00030101_H02_BBC3666 -- with the BBC3666 ID being in your PLINK data). The `idpat` specifies a regular expression that is used to extract out the ID you want. Two common patterns that are likely to be used are
   * `(.*)`  This is the default. This says that the ID as it appears in the sample sheet in the column "Institute Sample Label" _is_ the same as in your PLINK data. If you are so lucky, you need not do anything or even set `idpat` explicitly. 
   * `.*_(.*)` : this is for the case where the value in "Insitute Sample Label` is in the format PlateLabel_Well_SampleID and it says, ignore everything up to and including the last underscore and then everything else is the ID
   * If you have something else, you'll have to figure out how to either create a better sample sheet or pick the right regex. For Pythonites, the regular expression should contain exactly one pair of parenthesis which will be used by Python's regex features to extract out the ID.


### Performance parameters

There are three parameters that are important control performance. You probably won't need to change this, but feel free.

* `plink_mem_req` : specify in MB or GB how much memory your processes that use PLINK require;
* `other_mem_req` : specify how much other processes need;
* `max_cpus` : specify how many cores your PLINK processes can use. (This is only for those PLINK operations that are parallelisable. Some processes can't be parallelised and our workflow is designed so that for those processes only one core is used).


## 5. Output

A PDF report can be found in the output directory. This describes the process as well as what the inputs and outputs were.

Note that one issue that sometimes occurs in analysis is that there may over time be multple copies of the same file, perhaps with minor differences. To help version control, the PDF report captures the md5 checksums of inputs and outputs.


## 6. Example 
* data and command line can be found [h3agwas-examples](https://github.com/h3abionet/h3agwas-examples)

```
~/nextflow run h3abionet/h3agwas/main.nf --qc 1 --bfile data/array_plk/array --ouput_dir qc  --output kgpexample \
 --phenotype data/pheno/pheno_test.all --pheno_col phenoqc_ql \
 --case_control data/pheno/pheno_test.all --case_control_col Sex \
 --batch data/pheno/pheno_test.all --batch_col batch \
 -profile slurmSingularity -resume --qc  1
```

## 7. Overview
<img src="utils/qc_overview.png" title="overview of qc pipeline">


