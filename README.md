

<img src="aux/H3ABioNetlogo2.jpg"/>

# h3agwas Pipeline Version 2



## Background

h3aGWAS is a simple human GWAS analysis workflow for data quality control (QC) and basic association testing developed by [H3ABioNet](https://www.h3abionet.org/). It is an extension of the [witsGWAS pipeline](http://magosil86.github.io/witsGWAS/) for human genome-wide association studies built at the [Sydney Brenner Institute for Molecular Bioscience](https://www.wits.ac.za/research/sbimb/). h3aGWAS uses Nextflow as the basis for workflow managment and has been dockerised to facilitate portability.



The original version of the h3aGWAS was published in June 2017 with minor updates and bug fixes through the rest of the year. Based on experience with large data sets, the pipelines were considerably revised with additional features, reporting and a slightly different workflow.  


We have moved all scripts from Python 2 to Python 3, so you will need to have Python 3 installed.  

_Please ignore the Wiki in this version which refers to version 1_


## Brief introduction

A short video overview of the pipeline can be found at http://www.bioinf.wits.ac.za/h3a/h3agwas.mp4


### Restrictions
This version has been run on real data sets and works. However, not all cases have been thoroughly tested. In particular
* it is not robust when X chromosome data is not available
* the reporting assumes you want to do batch/site analysis. If you don't the code works but the report may look a bit odd with some figures repeated.
* we haven't tested fully with Singularity

The previous version 1 stable branch was commit bfd8c5a
(https://github.com/h3abionet/h3agwas/commit/bfd8c5a51ef85481e5590b8dfb3d46b5dd0cc77a)

There is one feature of the original workflow that has been omitted. Version 1 supported parallel GWAS anaysis of different data files in one Nextflow run. This has been removed. Although, not unuseful, this feature complicated the implementation and made expansion more difficulkt and also this capacity can be simulated easily at the operating system level.

The previous version has dependancies on Perl and R, which have been removed.

## Outline of documentation

1. Features
2. Installing the pipeline
3. A quick start example
4. The Nextflow configuration file
5. The QC pipeline: `plink-qc.nf`
6. A simple association testing pipeline: `plink-assoc.nf`
7. Converting Illumina genotyping reports to PLINK: `topbottom.nf`
8. Advanced options: Docker, PBS, Amazon EC2
9. Auxiliary Programs

# 1. Features

## 1.1 Goals of the h3agwas pipeline

The goals of this pipeline is to have a portable and robust pipeline
for performing a genome-wide association study

There are three separate workflows that make up *h3agwas*

1. `topbottom.nf`.  Conversion of Illumina genotyping reports with TOP/BOTTOM calls into PLINK format, aligning the calls.

2. `plink-qc.nf`: Quality control of the data. This is the focus of the pipeline. It takes as input PLINK data and has the following functions

   * Sample QC tasks checking:

       *  discordant sex information
       *  calculating missingness
       *  heterozygosity scores
       *  relatedness
       * discordant sex information
       * SNP QC tasks checking:

   * batch reports

       * remove duplicates
       * minor allele frequencies
       * SNP missingness
       * differential missingness
       * Hardy Weinberg Equilibrium deviations

3.  `plink-assoc.nf`: Association study. A simple analysis association study is done. The purpose of this is to give users an introduction to their data. Real studies, particularly those of the H3A consortium will have to handle compex co-variates and particular population study. We encourage users of our pipeline to submit thieir analysis for the use of other scientists.
  * Basic PLINK association tests, producing manhattan and qqplots
  * CMH association test - Association analysis, accounting for clusters
  * permutation testing
  * logistic regression
  * Efficient Mixed Model Association testing with gemma, boltlmm or fastlmm




## 1.2 Design principles

The goal of the H3ABionet GWAS pipeline is to provide a portable and robust pipeline for reproducble genome-wide association studies.


A GWAS requires a complex set of analyses with complex dependancies
between the analyses. We want to support GWAS work by supporting
* reproducibility -- we can rerun the entire analysis from start to finish;
* reusability -- we can run the entire analysis with different parameters in an efficient and consistent way;
* portability -- we can run the analysis on a laptop, on a server, on a cluster, in the cloud. The same workflow can be used for all environments, even if the time taken may change;


We took the following into account:

* The anticipated users are heterogeneous both in terms of their bioinformatics needs and the computing environments they will use.
* Each GWAS is different -- it must be customisable to allow the bioinformaticists to set different parameters.


There are two key technologies that we use, Nextflow and Docker, both of which support our design principles. 

### Nextflow

[Nextflow](https://www.nextflow.io/) is a workflow language designed at the Centre for Genomic Regulation, Barcelona. Although it is a general workflow language for science, it comes out of a bioinformmatics group and strongly supports bioinformatics. 

Our pipeline is built using Nextflow. However, users do not need to know anything about Nextflow. Obviously if you can do some programming you can customise and extend the pipelines, but you do not need to know Nextflow yourself. 

Nextlow is very easy to install and is highly portable. It supports partial execution and pipelines that scale.  Nextflow supports our worklow requirements very well.

### Docker

A GWAS requires several software tools to be installed. Using Docker we can simplify the installation. Essentially, Docker wraps up all software dependancies into _containers_. Instead of installing all the dependancies, you can install Docker, easily and then install our containers. (In fact you don't need to explicitly install our containers, Nextflow and our workflow will do that for you automatically).

We expect that many of our users will use Docker. However, we recognise that this won't be suitable for everyone because many high performance computing centres do not support Docker for security reasons. It is possible to run our pipeline without Docker and will give  instructions about which software needs to be installed.

### Singularity

The pipeline should support singularity. This is currently experimental though we hope to support it fully as a first-class citizen in the next release. We have run it through on a data set but we haven't tested all options and workflows.


# 2. Installing h3aGWAS

## 2.1 Background

The h3agwas pipeline can be run in different environments; the requirements differ. The different modes are described in detail below
* Running on Docker/Singularity. This is the easiest way of running h3agwas. We have a set of Docker containers that have all the required executables and libraries.
* Running natively on a local computer -- this is requires a number of external executables and libraries to be installed..
* Running with a scheduler -- Nextflow supports a range of schedulers. Our pipeline supports using docker or running natively.
* Running on Amazon EC2.  You need to have Amazon AWS credentials (and a credit card). Our EC2 pipeline uses Docker so this is very easy to run.
* We have also used Docker swarm. If you have a Docker swarm it's easy to do.

We now explore these in details

## 2.2 Basic requirements

**All** modes of h3agwas have the following requirements
* Java 8
* Nextflow. To install Nextflow, run the command below. It creates a _nextflow_ executable in the directory you ran the command. Move the executable to a directory on the system or user PATH and make it executable. You need to be running Nextflow 27 (January 2018) or later.
    `curl -fsSL get.nextflow.io | bash`
* Install the scripts
    `git clone https://github.com/h3agwas`

## 2.3 Installing with Docker

Install docker on your computer(s). Docker is available on most major platforms.  See [the Docker documentation](https://docs.docker.com/) for installation for your platform. 
    
That's it. 

## 2.4  Installing to run natively

This requires a standard Linux installation or macOS. It requires _bash_ to be available as the shell of the user running the pipeline.

The following code needs to be installed and placed in a directory on the user's PATH.

* plink 1.9 [Currently, it will not work on plink 2, though it is on our list of things to fix. It probably will work on plink 1.05 but just use plink 1.0]
* LaTeX. A standard installation of texlive should have all the packages you need. If you are installing a lightweight TeX version, you need the following pacakges which are part of texlive.: fancyhdr, datetime, geometry, graphicx, subfig, listings, longtable, array, booktabs, float, url.
* python 3.4 or later. pandas, numpy, matplotlib and openpyxl need to be installed. You can instally these by saying: `pip3 install pandas`  etc

If you want to run the `plink-assoc.nf` pipeline then you should install emmax and gemma if you are using those options.

# 3. Quick start example

This section shows a simple run of the `plink-qc.nf` pipeline that should run out of the box if you have installed the software or Docker. More details and general configuration will be shown later.

This section illustrates how to run the pipeline on a small sample data
file with default parameters.  For real runs, the data to be analysed
and the various parameters to be used are specified in the
_nextflow.config_ file.  The details will be explained in another
section.

The sample data to be used is in the _input_ directory (in PLINK
format as _sampleA.bed_, _sampleA.bim_, _sampleA.fam_). The default
_nextflow.config_ file uses this, and so you can run the workflow
through with this example. Note that this is a very small PLINK data set 
with no X-chromosome information and no sex checking is done.




## 3.1 Running on your local computer 

This requires that all software dependancies have been installed.

`nextflow run  plink-qc.nf`

The workflow runs and output goes to the _output_ directory. In the
_sampleA.pdf_ file, a record of the analysis can be found.

In order, to run the workflow on another PLINK data set, say _mydata.{bed,bim,fam}_, say

`nextflow run  plink-qc.nf --input_pat mydata`

If the data is another directory, and you want to the data to go elsehwere:

`nextflow run  plink-qc.nf --input_pat mydata --input_dir /data/project10/ --output_dir ~/results `

There are many other options that can be passed on the the command-line. Options can also be given in the _config_ file (explained below). We recommend putting options in the configuration file since these can be archived, which makes the workflow more portable

## 3.2 Running with Docker on your local computer

Execute 

`nextflow run  plink-qc.nf -profile docker`

Please note that the _first_ time you run the workflow using Docker,  the Docker images will be downloaded. *Warning:* This will take about 1GB of bandwidth which will consume bandwidth and will take time depending on your network connection. It is only the first time that the workflow runs that the image will be downloaded.


More options are shown later.


# 4 The Nextflow configuration file

Nextflow uses parameters that are passed to it and contents of a
configuration file to guide its behaviour. By default, the
configuration file used in _nextflow.config_. This includes
specifiying

* where the inputs come from and outputs go to;
* what the parameters of the various programs/steps. For example, in QC you can specify the what missingness cut-offs you want;
* the mode of operation -- for example, are you running it on a cluster? Using Docker?

To run your workflow, you need to modify the nextflow.config file, and
then run nexflow. Remember, that to make your workflow truly
reproducible you need to save a copy of the _config_ file. For this
reason although you can specify many parameters from the command line,
we recommend using the config file since this makes your runs
reproducible.  It may be useful to use git or similar tool to archive
your config files.

## 4.1 Specifiying an alternative configuration file

You can use the _-c_ option specify another configuration file in addition to the nextflow.config file

```nextflow run -c data1.config plink-qc.nf```


**This is highly recommended.** We recommend that you keep the `nextflow.config` file as static as possible, perhaps not even modifying it from the default config. Then  for any
 run or data set, have a much smaller config file that only specifies the changes you want made. The base `nextflow.config` file will typically contain config options that are best set by the h3aGWAS developers (e.g., the names of the docker containers) or default GWAS options that are unlikely to change. In your separate config file, you will specify the run-specific options, such as data sets, directories or particular GWAS parameters you want. Both configuration files should be specified. For example, suppose I create a sub-directory within the directory where the nextflow file is (probably called h3agwas). Within the h3agwas directory I keep my nexflow.config file and the nextflow file itself. From the sub-directory, I run the workflow by saying:

```nextflow run  -c data1.config ../plink-qc.nf```

This will automatically use the `nextflow.config` file in either the current or parent directory. Note that the the config files are processed in order: if an option is set into two config files, the latter one takes precedence.


## 4.2 Creating an auxiliary nextflow .config file

There is a template of a nextflow.config file called aux.config.template. This is a read only file. Make a copy of it, call it _aux.config_ (or some suitable name).  This file contains all the options a user is likely to want to change. It does not specify options like the names of docker containers etc. Of course, you can if you wish modify the nextflow.config file, but we recommend against it. Your auxiliary file should supplement the nextflow.config file.

Then fill in the details in the config that are required for your run. These are expained in more detail below.

## 4.3 Using the Excel spreadsheet template

For many users it may be convenient to use the Excel spreadsheet (config.xlsx and a read-only template file config.xlsx.template). This can be used just as an _aide-memoire_, but we also have an auxiliary program that converts the Excel spreadsheet into a config file. The program _config-gen/dist/config-gen.jar_ takes the spreadsheet and produces a config file.

The spreadsheet has the following columns
* A. a brief one-line description of the parameter;
* B. the name of the parameter as found in the config file;
* C. the default value that will be used by  the _config-gen_ program if no value specified in column E;
* D. possibie alternate value the user might consider; 
* E. the value that the user wants to use

If you are using this semi-automated way of producing the config file, remember that to be fully reproducible the config file must be saved too. We suggest making a copy of the spreadsheet template file giving it an appropriate name.

You need to enter something into column E.

Remember this creates only the auxiliary config file -- you still need the main file

To run the _config-gen_ program, you would for example say:

`java -jar ./config-gen/dist/config-gen.jar nameofspreadsheet > run10.config`

I suggest editing that by removing any lines you don't want to change.

Then you would run your script by saying

`nextflow run -c run10.config plink-qc.nf`

The _nextflow.config_ file will automatically be used, except for any additions or changes that are in the _run10.config_ file.

## 4.4 Specifying options

When you run the the scripts there are a number of different options that you might want to use. These options are specified by using  the `-flag` or `--flag` notation. The flags with a single hyphen (e.g. `-resume`) are standard Nextflow options applicable to all Nextflow scripts. The flags with a double hyphen (e.g., `--pi_hat`) are options that are specific to _our_ scripts.  *Take care not to mix this up as it's an easy error to make, and may cause silent errors to occur.*


Almost all the workflow options that are in the _nextflow.config_ file can also be passed on the command line and they will then override anything in the config like. For example

```nextflow run plink-qc.nf   --cut_miss  0.04```

sets the maximim allowable per-SNP misisng to 4%. However, this should only be used when debugging and playing round. Rather, keep the options in the auxiliary config file that you save. By putting options on the command line you reduce reproducibility. (Using the parameters that change the mode of the running -- e.g. whether using docker or whether to produce a time line only affects time taken and auxiliary data rather than the substantive results).


## 4.5 Partial execution and resuming execution

Often a workflow may fail in the middle of execution because there's a problem with data (perhaps a typo in the name of a file), or you may want to run the workflow with slightly different parameters. Nextflow is very good in detecting what parts of the workflow need to re-executed -- use the `-resume` option. 

## 4.6 Cleaning up 

If you want to clean up your work directory, say `nextflow clean`.

## 4.7 Workflow overview, and timing

Nextflow provides [several options](https://www.nextflow.io/docs/latest/tracing.html) for visualising and tracing workflow. See the Nextflow documentation for details. Two of the options are:

* A nice graphic of a run of your workflow

    `nextflow run plink-qc.nf -with-dag quality-d.pdf`

* A timeline of your workflow and individual processes (produced as an html file).

    `nextflow run <pipeline name> -with-timeline time.html`

    This is useful for seeing how long different parts of your process took. Also useful is peak virtual memory used, which you may need to know if running on very large data to ensure you have a big enough machine and specify the right parmeters.





# 5 The QC pipeline: `plink-qc.nf`


This section describes the various ways in which the pipeline can be run and various options. Usually options are specified in the _nextflow.config_ file (or which ever file you use). However, you can also pass parameters to the Nextflow script on the command-line. Parameters on the command line over-ride any parameters specified in the config file.


The main pipeline is the PLINK QC pipeline. It takes as input PLINK bed,bim,fam files and performs quality control on  the data according to the parameters specified in the config file.

The Nextflow script file is called *plink-qc.nf*. This could be
called, for example, by running `nextflow run plink-qc.nf`.

The output of the QC is a set of PLINK files that can be used for GWAS, as well as PDF report that describes the QC steps.

## 5.1 Input/Output :  PLINK format

Users will run the pipeline giving as input PLINK 1.9 bed, bim and fam files.  The key Nextflow parameters to set are:

* `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
* input, output and script directories: the default is that these are subdirectories of the `work_dir` and there'll seldom be reason to change these;
* `input_pat` : this typically will be the base name of the PLINK files you want to process (i.e., do not include the file suffix). But you could be put any Unix-style glob here. The workflow will match files in the relevant `input_dir` directory;
* `high_ld_regions_fname`: this is optional -- it is a list of regions which are in very high LD -- and are exclude when checking for relationships (https://www.cog-genomics.org/plink/1.9/filter#mrange_id).  Provide either absolute file path or relative to where you are running. In a previous version this was relative to input_dir, which is not right.
See [https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)](https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)) for a discussion.

* `output`: the base name of the output files. *This cannot be the same as the input!!!*

## 5.2 Overview of the workflow

The QC process consists of:

* removing duplicate markers;
* indentifying indviduals for whom there is discordant sex information;
* removing individuals with too high missingness or excessive heterozygosity;
* detecting whether there are any related individuals and removing enough to ensure that there are not related pairs;
* removing SNPs with too low MAF, or too high missingness, or anomalous HWE, or SNPs where there is a high differential missingness between cases and controls;
* a PCA of the resultant data is computed;
* a detailed report of the QC process is done.

## 5.3 Additional QC Parameters

The following parameters control QC

*  `sexinfo_available`: `true` or `false`. If we don't have sex information then we cannot do the check for discordant genotype. Note that it does not make sense (and is an error) to have sexinfo_available set to true when there is no X-chromosme data in the file;
*  `f_low_male` and `f_hi_female`. Discordant sex genotype is done on the X-chromosome using the non-recombining parts. F, the in-breeding coefficient of the X-chromosome is computed. If F is above `f_low_male`, the individual is predicted to be male; if F is below `f_hi_female`, the individual is predicted to be female. Anyone in between is flagged. These cut-off values are arbitrary and especially in large samples you are likely to find a range of F values. However, a large number of discrepant values probably indicates a sample mishandle error.  The PLINK default values (0.8 and 0.2) are the default parameters of the pipeline.
*  `cut_het_high`: What is the maximum allowable heterozygosity for individualsl;
*  `cut_het_low`: minimum
*   `cut_maf `: the minimum minor allele frequency a SNP must have to be included
*   `cut_diff_miss `: allowable differential missingness between cases and controls;
*   `cut_geno`: maximum allowable per-SNP mssingness
*   `cut_mind`: maximum allowable per-individual missingness
*   `cut_hwe`: minimum allowable per-SNP Hardy-Weinberg Equilibrium p-value 
*   `pi_hat`:  maximum allowable relatedness
*   `remove_on_bp`: the first step in the pipeline is to remove duplicate SNPs. There are two ways of detecting duplicates. First, if SNPs have duplicate names (column 1 -- numbering from 0 -- of the bim file). We always remove duplicate SNPs based on this since PLINK gets very upset otherwise. Second, if they are at the same chromosome and base position. If this variable is set to 1, then duplicates based on chromosome or base position are removed too.
*   `batch`: if you want to do QC at a batch level, you need to specify a file with the batch information. This should be a standard PLINK phenotype file (with labels as the first line). If you specify "false" or 0, then no batch-analysis is done. Typically batch information relates to the batches in which samples were genotyped is not intrinsic to the data (e.g. you genotype the first 2500 samples that are available).
*   `batch_col`: the column label of the file to be used.
*   `phenotype`: default is 'false'. If you are doing batch analysis you may wish to show how different sub-groups perform in QC with respect to the batch. You will then specify a PLINK-style phenotype file (with labels as the first name).  For example, if you have a multi-site project, you may choose to use the site information as a phenotype. Other possibilities are sex and self-identified group. If you specify "false" or 0, no categorisation will be done.
* `pheno_col` is the column label of the column in  the phenotype file which should be used.
*  `case_control` : This is the name of a PLINK-style phenotype file with labels in the first line. This is a compulsory parameter. The QC process uses the case/control status of individuals. A principal component analysis is done. We do not expect typically overall that there will be difference between cases and controls. The PC-analysis tests that this is so. Of course, you need to apply your mind to the result as YMMV. If your study has several case/control categories, choose an appropriate one that will give insight. If you only have continuous measures (e.g., BMI), then discretise and make an artificial case-control category. Remember, this is for QC purposes not to find interesting biology.
* `case_control_col`: this is the label of the column.

Several of the above parameters make reference to a phenotype file. Of course, these can be to the same phenotype file, but probably using different columns.

## 5.4 Performance parameters

There are three parameters that are important control performance. You probably won't need to change this, but feel free.

* `plink_process_memory` : specify in MB or GB how much memory your processes that use PLINK require;
*  ` other_process_memory` : specify how much other processes need;
*  `max_plink_cores` : specify how many cores your PLINK processes can use. (This is only for those PLINK operations that are parallelisable. Some processes can't be parallelised and our workflow is designed so that for those processes only one core is used).


## 5.5 Output

A PDF report can be found in the output directory. This describes the process as well as what the inputs and outputs were.

Note that one issue that sometimes occurs in analysis is that there may over time be multple copies of the same file, perhaps with minor differences. To help version control, the PDF report captures the md5 checksums of inputs and outputs.


# 6. Simple association test pipeline: `plink-assoc.nf`

An association study is a complex analysis and each analysis has to consider
* the disease/phenotype being studied and its mode of inheritance
* population structure
* other covariates

For this reason it is difficult to build a high quality, generic pipeline to do an association study. 

The purpose of this pipeline is to perform a very superficial initial analysis that can be used as one piece of information to guide a rigorous analysis. Of course, we would encourage users to build their own Nextflow script for their rigorous analysis, perhaps using our script as a start.

Our script, *plink-assoc.nf* takes as input PLINK files that have been through quality control and 
* does a principal component analysis on the data, and produces pictures from that; 
* performs a simple association test giving odds ratio and  raw and adjusted _p_ values

## Running

The pipeline is run: `nextflow run plink-assoc.nf`

The key options are:
* `input_dir`, `output_dir`: where input and output goes to and comes from;
* `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
* `data`: a tab-separated file that contains phenotype data for your particpants. The row is a header line, with one participant per line after that. The first two columns should FID IID followed by any phenotype values you want to use for testing or for covariates. It can contain other data too -- as long as the ones that you need are in this file.
* `pheno`: a comma-separated list of phenotypes that you want to test. Each phenotype is tested separately. If you wish to do a transform to the phenottype, you can suffix the phenotype with _/function_ where function is a numpy function. For example you might say `--pheno bmi/np.log` which will apply the log function to the phenotype. Any numpy function can be used, typical uses being np.log and np.sqrt. We plan to support user provision of a user-given function.
* `covariates`: a comma-separated list of phenotypes that you want to use
* `exclude_snps` option to exclude some snps active just for boltlmm (see `--exclude` in boltlmm manual) : SNP ID must be first column (default none)
*  `print_pca` : by default pipeline compute and print pca (`print_pca`=1), if you want avoid this step (`print_pca` = 0)
*  `file_rs_buildrelat` : file with rs list (one by lines) to build genetics models (relatdness), for gemma `-snps` for boltlmm `--modelSnps`
* `genetic_map_file` : genetic map used in boltlmm 

By default a chi2 test for association is done. But you can do multiple different tests in one run by settintg the appropriate parameter to 1. Note at least one must be set to 1

 * `chi2` : should a chi2 test be used (0 or 1)
 * `fisher`: Fisher exact test
 *  `linear`: should linear regreession be used?
 *  `logistic`: should linear regression be used?
 *  `gemma`: should gemma be used?
    *  see [manual](www.xzlab.org/software/GEMMAmanual.pdf)
    *  `gemma_num_cores`: if gemma is used set this up to 8
    *  `gemma_mem_req`: For 10k samples, 2 million SNPs, we needed 4GB of RAM (default : 6GB)
    *  `gemma_mat_rel` : file contains in gemma format matrix relatdness used by gemma  (format 1, see manual), matrix must be in same order than fam file. Avoid computation of relatdness by pipeline. 
 *  `boltlmm`: should boltlmm be used? 
    *  see [manual](https://data.broadinstitute.org/alkesgroup/BOLT-LMM/)
    * if SNPs is higher than 950,000, 950,000 SNPs are chosen randomly to build the model (see --modelSnps option in bolt)
    * `bolt_covariates_type` : for bolt need to define if covariate is binary (0) or continue (1), a comma-separated list as same order than covariates 
    * `bolt_ld_score_file` : A table of reference LD scores for boltlmm is needed to calibrate the BOLT-LMM statistic (option in boltlmm --LDscoresFile),to choice a specific column in Ld file you can use `bolt_ld_scores_col` option (by default : LDSCORE) if option is not provided --LDscoresUseChip used.
    * `bolt_use_missing_cov` : option to "missing indicator method", by default no activate (0), to activate (1) (--covarUseMissingIndic option in boltlmm), which adds indicator variables demarcating missing status as additional covariates.
    * `bolt_num_cores` if bolt is used set this up to 8
    * `bolt_mem_req` memory required for boltlmm, (default : 6GB)
    * impute2 data in bolt  :
      * bolt_impute2filelist : list of impute2 files, each line contains : `chronumber` `file`, file must be in full pattern
      *`bolt_impute2fidiid` : list of individual in same order than bolt_impute2filelist
 *  `fastlmm`: should fastlmm be used?
    *  see [manual](https://github.com/MicrosoftGenomics/FaST-LMM)
    * `fastlmm_num_cores`: if fastmll is used set this up to 8
    * `fastlmm_mem_req`: memory required for fasttlmm (default : 15GB)
    * `fastlmm_multi` : memory used by fastmll is very big and increase with snp number, option run fastlmm by chromosome, with relatedness matrix computed before with gemma (-gk 1) 
    * `fastlmmc_bin` : should change binary for fastlmmc (default fastlmmc)
and then for all the tests except _gemma_, _boltlmm_ and _fastlmm_, do you want to adjust for multiple testing 

* `adjust`: do we want to do explicit testing for Bonferroni correction et al that PLINK odes
* `mperm`: do you want to test doing permutation testing. If so, how many tests?  By default this is 1000.

with pipeline, do a GxE interaction with Gemma and Plink, arguments :
  * `gxe` : environmental variables to do gxe analysis with `pheno`, must be coded in 1 and 2 for plink
  * `gemma_gxe` : GxE interation with gemma [default : 0], see  `covariates` to add covariates in gemma models
  * `plink_gxe` : GxE interation with plink (see option -gxe, in [plink manual](http://zzz.bwh.harvard.edu/plink/anal.shtml#qtgxe)) [default : 0], no covariate could be provided.


For example

```nextflow run plink-assoc --input_pat raw-GWA-data --chi2 1 --logistic 1 --adjust 1```

analyses the files `raw-GWA-data` bed, bim, fam files and performs a chi2 and logistic regression test, and also does multiple testing correction.

Other flags are:
* `thin`. You can set this to a floating point number in the range (0, 1] and then the PLINK data files are thinned leaving only that proportion of the SNPs. This allows pipeline to be tested with a small proportion of the data This is probably only needed for debugging purposes and usually this should not be be set.
* `chrom`. Only do testing on this chromosome.

# 7. Converting from Illumina genotyping reports in TOP/BOTTOM format

This workflow is run by 

```nextflow run topbottom.nf```

and converts from an Illumina TOP/BOTTOM call file. Together with auxiliary input data, this file is first converted into a raw PLINK file and then the PLINK file is aligned to a strand, and then convered into binary PLINK format. This process can take a very long time.


This process is expensive because:

* the top/bottom file is a very bulky and inefficient format
* we convert first to PLINK using the inefficenct PED format

As an example, on a 2.5m SNP-chip with 10k individuals, you are looking at over 200 CPU-hours.

17 January 2017: this code has been completely reworked to make it more efficient and with fewer dependancies. But it is also less  powerful. There is other code, such as Don Armstrong's code and another option is the unofficially supported Illumina code https://github.com/Illumina/GTCtoVCF.  You can go back to the older code

## Input 

You require the following input
* the actual call files from Illumina
* `input_dir`: the directories where the Illumina genotyping reports can be found. Unix-style globs allowed

e.g. `params.input_dir = "/project/HumCVD/Batches/Batch*/Batch*Reports/"`

* `input_pat`. The files that inside these directories.

e.g. `params.input_pat = "*gtReport*.csv.gz"`

* `output`: the base name of the PLINK output file

e.g. `params.output = "cvd-rawcalls"`

* `manifest`: The chip manifest file. This is crucial. You can find examples here: https://support.illumina.com/downloads.html, but you may have to ask Illumina for the correct vrsion.

* `chipdescription`: this is a standard file produced by Illumina for your chip which contains not only the chromosome/coordinate of each SNP but also the genomic position (measured in centimorgans). If you don't have this -- give the manifest file. All will work except your bim files will not contain genonomic positoin

* `samplesheet`: This is Excel spreadsheet that Illumina provides which details each perfson in the study for whom you have genotyping results. If you don't have it, you can set this variable to 0 or the empty string, in which case the output PLINK fam file will have unknown values for sex and phenotype.  Alternatively, ifyou don't have it, you can make your own. There are three columns that are important: "Institute Sample Label", "Manifest Gender" and "Batch Comment". These must be there. The _label_ is the ID of the person. In the current workflow this ID is used for both the FID and IID. If you have a family study you may need to manually change the fam file.

Please note that *we expect all entries in the sample IDs etc to be alphanumeric 0-9, Latin letters (NO accents!), underscore, space, hyphen*. The code may break otherwise.

* `idpat`. Default is "0" or "" (ignore). By default, we use the sample ID as found in the genotype report and sample sheet. PLINK fam files require a double barrelled name (FID IID) -- we just double the ID as found. However, this may not be ideal since
the Illumina IDs in the sample ID are typically a long string some of  the components of which will not be useful when you are analysing the result. You can change the sample ID by providing a Python-style regular expression which decribes the components. The regex groups describe the components. If there is one group, it is doubled. If there are two groups, then those become the FID and IID. Only one or two groups are permissible. 

For example, suppose the ID as found in the Illumina input data is `WG0680781-DNA_A02_ABCDE`, if you use ".*_(.+)" as the idpat, then the FID IID used would be ABCDE ABCDE. If you used "(\\w+)_DNA_(\\w+)_" then the FID IIS used would be "WG0680781 A02". Note how we need to escape the backslash twice.

* `output_align`. This can be one of three values: _topbottom_, _dbsnp_ and _ref_. If topbot, the SNPs will be aligned to the Illumnia TOP strand. If dbsnp, the output will be aligned to the dbSNP report, if "ref", the output will be aligned to a given reference strand. In the latter two cases, many of the SNPs will be flipped (e.g. an A/C SNP will become G/T; and A/T SNP will become T/A).

* `strandreport`: This is an Illumina-style strand report. It is not needed if you choose "tobot" above, but it is needed for "dbsnp" or "ref"

* `refererence`: This is the name of a file that gives the reference allele for each SNP on the chip.  This is only useful if the "ref" option is used for `output_align`, and is optional in this case. Note that the difference between aligning and the use of the reference. Aligning will decide which strand of the reference genome as found in the Illumina genotyping teh alleles can be found on. For example, if the genotyping report gives the two options as A and C, aligning checks whether this is A and C on the + strand of the reference genome or (and so will be A and C in the output bim file) or whther this is A and C on the $-$ strand of the reference genome and so should be reported as T and G. This is done using information in the chip manifest file. The second step is to know which allele is the reference allele and which is the alternate allele.

A reference file suitable for the H3A chip can be found here http://www.bioinf.wits.ac.za/data/h3agwas/. Two formats for the reference file are supported: (1) simple -- two columns, no header row, the first column ins the SNP ID, the second the reference allele; and (2) complex -- >= two columns, one header row, the header row must contain a column label SNP and a column label Base for the SNP ID and reference allele respectively, all other columns are ignored.

* `batch_col`: For this workflow, the `batch_col` parameter is a column in the `samplesheet` that should be used to extract out out the 6-th column of the `fam` file, or the phenotype. This allows you do do batch analysis. Of course, you can choose anything you like to be the "batch". The default value is 0, which means just set the 6-th column of the fam file to -9.  One special case: If the contents of the column is of the form "Batch n", then only the _n_ is returned.

* `samplesize`: This was  included mainly for development purposes but _perhaps_ might be helpful to some users. This allows you sample only the first _n_ people in each genotype report. This allows you to extract out a small subset of the data for testing purposes. The default is 0, which means that *all* individuals will be generated.

## Output

The output are a set of PLINK files (bed, bim, fam, log).

In addition, there may be a file with a _.badsnps_. If you chose to align against the reference genome, these are the SNPs for which the reference allele is inconsistent with the two allele choices in the data. For example, the reference allele is A and the choice of alleles in the data is C/T.  Hopefully this will be a small number (a thousand or so). There are a number of reasons why this may be the case;

* There is a problem in the chip. Chip design and some SNPs are hard to design for. This may result in some SNPs coordinates that mismatch with the probe design, or for which the strand alignment is unclear.
* There is a problem with the reference file you provided.

Possible ways forward:
* if you have tens of thousands of such SNPs, then there's a problem which must be resolved.
* if you have one or two thousand or less, then for population structure studies, just remove the SNPs from the data.
* if you are doing a GWAS, you can leave them in but if you find interesting matches check carefully whether the SNP is on this list. If it is you need to be carefully study the SNP (looking at the reference genome, strand file, etc to see why it mismatched and looking at the image files).
* if you are doing imputation from the data, it's probably best to remove them unless you have a specific need for a particular SNP.


# 8. Running the workflow in different environments

In the  quick start we gave an overview of running our workflows in different environments. Here we go through all the options, in a little more detail

## 8.1 Running natively on a machine

This option requires that all dependancies have been installed. You run the code by saying

```
nextflow run plink-qc.nf
```

You can add that any extra parameters at the end.

## 8.2 Running  on a local machine with Docker

This requires the user to have docker installed.

Run by `nextlow run plink-qc.nf -profile docker`



## 8.3 Running on a cluster 

Nextflow supports execution on clusters using standard resource managers, including Torque/PBS, SLURM and SGE. Log on to the head node of the cluster, and execute the workflow as shown below. Nextflow submits the jobs to the cluster on your behalf, taking care of any dependancies. If your job is likely to run for a long time because you've got really large data sets, use a tool like _screen_ to allow you to control your session without timing out.

To run using Torque/PBS, log into the head node. Edit the _nextflow.config_ file (either directly or using our helper script). If you are doing this manually, edit the _nextflow.config_ file by looking for the stanza (around like 85)


```
    pbs {
        process.executor = 'pbs'
	process.queue = 'long'
    }
```
and change _long_ to whatever queue you are using. Note that in the current version, the only way in which the queue can be changed is by manuall editing the nextflow.config file. You can either change the `process.queue` line, or, better, modify the definition of the queue variable at the top of the `nextflow.config` file.


Then you can run by saying

```
nextflow run plink-qc.nf -profile pbs
```

If you are using another scheduler, the changes should be straight-forward. For example, to run using SLURM, add a stanza like within the _profile_environment of the _nextflow.config_ file 

```
    slurm {
        process.executor = 'slurm'
	process.queue = 'long'
    }
```

and  then use this as the profile.

## Running on a cluster with Docker

If you have a cluster which runs Docker, you can get the best of both worlds by editing the queue variable in the _pbsDocker_ stanza, and then running

```
nextflow run plink-qc.nf -profile pbsDocker
```

We assume all the data is visible to all nodes in the swarm. Log into the head node of the Swarm and run your chosed workflow -- for example


## 8.4 Running on Docker Swarm

We have tested our workflow on different Docker Swarms. How to set up Docker Swarm is beyond the scope of this tutorial, but if you have a Docker Swarm, it is easy to run. From the head node of your Docker swarm, run

```
nextflow run plink-qc.nf -profile dockerSwarm
```

## 8.5 Singularity

The workflows run on Singularity thought this is currently experimental and we haven't  tried and tested all options. We don't have first class support for Singularity yet, so you will have do some PT but it's not too bad.

### Get the Singularity images

You need to make get the Singularity images for the workflow you want. There are two options

* running `docker2singularity` (https://github.com/singularityware/docker2singularity)

This is a good option for people who have a computer with docker running (e.g., their desktop) and will then move the Singularity container to another computer which doesn't run cluster.

An example run would be to create a directory _singularity_ somewhere and then run the _docker2singularity_ workflow. The singularity image will be put in the specified directory. For example, I did:

```
  docker run -v /var/run/docker.sock:/var/run/docker.sock  \  
           -v  /Users/scott/singularity:/output --privileged -t --rm  \
           singularityware/docker2singularity quay.io/h3abionet_org/py3plink 
  docker run -v /var/run/docker.sock:/var/run/docker.sock  \  
             -v  /Users/scott/singularity:/output --privileged -t --rm  \
             singularityware/docker2singularity quay.io/h3abionet_org/h3agwas-texlive
```

* Using `singularity pull`

This is easier, but the images are bigger

```
singularity pull --size 1880  docker://quay.io/h3abionet_org/py3plink 
singularity pull --size 1880  docker://quay.io/h3abionet_org/h3agwas-texlive

```

### Copy the images

Move the images to the system you want to run the workflow on. If you're on a cluster, then this must be on a system-wide  file system

### Edit the nextflow.config file

You need to edit this part of the _singularity_ stanza

```
        sg_py3Image = "/home/scott/py3plink.img"
        sg_latexImage = "/home/scott/h3agwas-texlive.img"

        process.executor = 'pbs'
        process.queue = 'batch'
```

The two image variables should be set to where you have put your singularity images. The `process.executor` variabel should be set to `local` if you want to run the workflow on the local computer, and to `pbs` if on a cluster using PBS. In the latter case, you should also set the queue variable appropriatelly.


## Run the workflow

`nextflow run plink-qc.nf --profile singularity`


## 8.5 Other container services

We hope to support Singularity soon. We are unlikely to support udocker unless Nextflow does. See this link for a discussion https://www.nextflow.io/blog/2016/more-fun-containers-hpc.html

## 8.6 Running on Amazon EC2


Nextflow supports execution on Amazon EC2. Of course, you can do your own custom thing on Amazon EC2, but there is direct support from Nextflow and  we provide an Amazon AMI that allows you to use Amazon very easilyl. This discussion assumes you are familiar with Amazon and EC2 and shows you how to run the workflow on EC2:


1. We assume you have an Amazon AWS account and have some familiariy with EC2. The easiest way to run is by building am Amazon Elastic File System (EFS) which persists between runs. Each time you run,  you attach the EFS to the cluster you use. We assume you have
* Your Amazon accessKey and secretKey
* you have the ID of your EFS
* you have the ID of the subnet you will use for your Amazon EC2.

    Edit the  nextflow config file to add your keys to the _aws_ stanza, as well as changing the AMI ID, sharedStorageID, the mount and subnet ID. *BUT see point 8 below for a better way of doing things*.
    ```
    aws {
       accessKey ='AAAAAAAAAAAAAAAAA'
       secretKey = 'raghdkGAHHGH13hg3hGAH18382GAJHAJHG11'
       region    ='eu-west-1'
    }

    cloud {
             ...
             ...
             ... other options
             imageId = "ami-710b9108"      // AMI which has cloud-init installed
             sharedStorageId   = "fs-XXXXXXXXX"   // Set a common mount point for images
             sharedStorageMount = "/mnt/shared 
   	     subnetId = "subnet-XXXXXXX" 
    }
     ```

    Note that the AMI is the H3ABionet AMI ID, which you should use. The other information such as the keys, sharedStorageID and subnetID you have to set to what you have.

The instructions below assume you are using nextflow. If you launch the machine directly, the user will be `ec2-user`; if you use the instructions below, you will be told who the user on Amazon instance is (probably the same userid as your own machine).

2. Create the cloud. For the simple example, you only need to have one machine. If you have many, big files adjust accordingly.

    `nextflow cloud create h3agwascloud -c 1`
 
    The name of the cluster is your choice (_h3agwascloud_ is your choice).
    

3.  If successful, you will be given the ID of the headnode of the cluster to log in. You should see a message like,

```
> cluster name: h3agwascloud
> instances count: 1
> Launch configuration:
 - bootStorageSize: '20GB'
 - driver: 'aws'
 - imageId: 'ami-710b9108'
 - instanceType: 'm4.xlarge'
 - keyFile: /home/user/.ssh/id_rsa.pub
 - sharedStorageId: 'fs-e17f461c'
 - sharedStorageMount: '/mnt/shared'
 - subnetId: 'subnet-b321c8c2'
 - userName: 'scott'
 - autoscale:
   - enabled: true
   - maxInstances: 5
   - terminateWhenIdle: true



Please confirm you really want to launch the cluster with above configuration [y/n] y
Launching master node -- Waiting for `running` status.. ready.
Login in the master node using the following command: 
  ssh -i /home/scott/.ssh/id_rsa scott@ec2-54-246-155-85.eu-west-1.compute.amazonaws.com
```

4. ssh into the head node of your Amazon cluster. The EFS is mounted onto `/mnt/shared`. In our example, we will analyse the files _sampleA.{bed,bim,fam}_ in the /mnt/shared/input directory  The  _nextflow_ binary will be found in your home directory. (Note that you can choose to mount the EFS on another mount point by modifying the nextflow option `sharedStorageMount`;

5. For real runs, upload any data you need. I suggest you put in the /mnt/shared directory, and do not put any data output on the home directory. Yo

5. Run the workflow -- you can run directly from github. The AMI doesn't have any of the bioinformatics software installed. 

    Specify the docker profile and nextflow will run using Docker, fetching any necessary images.

    Do `nextflow pull h3abionet/h3agwas`

    This pull is not strictly necessary the first time you run the job, but it's a  good practice to get into to check if there are updates.

6. Then run the workflow

    `nextflow run  h3abionet/h3agwas -profile docker --input_dir=/mnt/shared/XXXXX/projects/h3abionet/h3agwas/input/ --work_dir=/mnt/shared `

   You will need to replace XXXXX with your userid -- the local copy of the repo is found in the `/mnt/shared/XXXXX/projects/h3abionet/h3agwas/` directory. But we want the work directory to be elsewhere.

   Of course, you can also use other parameters (e.g. -resume or --work_dir). For your own run you will want to use your nextflow.config file.


   By default, running the workflow like this runs the `plink-qc.nf` script. If you want to run one of the other scripts you would say `nextflow run  h3abionet/h3agwas/topbottom.nf` or `nextflow run h3abionet/h3agwas/plink-assoc.nf` etc. 


7. The output of the default runcan be found in` /mnt/shared/output`. The file sampleA.pdf is a report of the analysis that was done.

8. Remember to shutdown the Amazon cluster to avoid unduly boosting Amazon's share price.

    `nextflow cloud shutdown h3agwascloud`


9. _Security considerations_: Note that your Amazon credentials should be kept confidential. Practically this means adding the credentials to your _nextflow.config_ file is a bad idea, especially if you put that under git control or if you share your nextflow scripts. So a better way of handling this is to  put confidential information in a separate file that you don't share. So I have a file called _scott.aws_
which has the following:
```
aws {
    accessKey ='APT3YGD76GNbOP1HSTYU4'
    secretKey = 'WHATEVERYOURSECRETKEYISGOESHERE'
    region    ='eu-west-1'
}

cloud {
            sharedStorageId   = "fs-XXXXXX"   
	    subnetId = "subnet-XXXXXX" 
}

```
Then when you create your cloud you say this on your local machine
  
  `nextflow -c scott.aws -c run10.config cloud create scottcluster -c 5`

Note there are two uses of `-c`. The positions of these arguments are crucial. The first is an argument to _nextflow_ itself and gives a configuration file to nextflow to use. The second is an argument to _cloud create_ which says how many nodes should be created. 

The _scott.aws_ file is not shared or put under git control. The _nextflow.config_ and _run10.config_ files can be archived, put under git control and so on because you _want_ to share and archive this information with o thers.




# 9. Auxiliary Programs

These are in the aux directory

## 9.1 updateFam.py

Can be used to update fam files. You probably won't need it, but others might find it useful. The intended application might be that there's been a mix-up of sample IDs and you want to correct.  The program takes four parameters: the original sample sheet, a new sample sheet (only has to include those elements that have changed), the original fam file, and then the base of a newfam file name.  The program takes the plate and well as the authorative ID of a sample. For every row in the updated sheet, the program finds the plate and well, looks up the corresponded entry in the original sheet, and then replaces that associated ID in the fam file. For example, if we have

_Original sheet_
```
Plate   Well Sample 
W77888  G01  AAAAAA
```

_New sheet_
```
Plate   Well Sample 
W77888  G01  BBBBBB
```

Then the new fam file has the AAAAA entry replaced with the BBBBB entry

Three files are output: a fam file, an error file (the IDs of individuals who are in th e sample sheet but not the fam file are output), and a switch file (containing all the changes that were made). Some problems like duplicate entries are detected.

## 9.2 getRunsTimes.pl (By Harry Noyes)

Nextflow has great options for showing resourc usage. However, you have to remember to set those option when you run.  It's easy to forget to do this. This very useful script by Harry Noyes (harry@liverpool.ac.uk) parses the .nextflow.log file  for you

# 10 Simulation pipeline: `simul-assoc.nf`

This section describes a pipeline in devlopment, purpose of this pipeline is to estimate false positive and false negative with simulated phenotype, Our script, *simul-assoc.nf* takes as input PLINK files that have been through quality control and
  * Simulate quantitative phenotypes with [phenosim](https://www.ncbi.nlm.nih.gov/pubmed/21714868) based on genetics data 
  * perform a GWAS on  phenotype simulated using gemma, boltlmm.
  * Perform summary statistics.

## Installation
a version of _phenosim_ adapted is already in nextflow binary, write in python2. plink, gemma and bolt must be installed 

## Running

The pipeline is run: `nextflow run simul-assoc.nf`

The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * input, output and script directories: the default is that these are subdirectories of the `work_dir` and there'll seldom be reason to change these;
  * `input_pat` : this typically will be the base name of the PLINK files you want to process (i.e., do not include the file suffix). But you could be put any Unix-style glob here. The workflow will match files in the relevant `input_dir` directory;
  * num_cores : cores number used 
  * ph_mem_req : memory request for phenosim
  *  Simulation option :
     * `phs_nb_sim` : simulation number (default : 5) 
     * `phs_quant_trait` :  quantitative trait simulation : 1, qualitative not develop yet (default : 1, -q option in phenosim)
     * Quantitative trait option :
        * `ph_nb_qtl` : number of simulated QTN (default: 2, option -n in phenosim)
        * `ph_list_qtl` : proportion of variance explained by each QTNs, separate the values by commas (default : 0.05 -q in phenosim)
        * `ph_maf_r` :  MAF range for causal markers (upper and lower bound, separated by a comma, no space) (default: 0.05,1.0, -maf_r in phenosim)
        * option to do a linear transformation of phenotype with co factor of external data and normatisation:
           * ph_normalise : perform a normalisation (1) or not 0 (Default)
           * each phenotype i be normalise using newpheno = norm(pheno)+var0i*a+var1i*b+ ... + intercept
           * `ph_cov_norm` : contains coefficients for relation separed by a comma (ex "sex=0.2,age=-0.1)
           * `data` : contains cofactor data for each individuals used to normalise with 
           * `ph_cov_range` : normalisation range for initial phenotype
           * `ph_intercept` : intercept
  * Association option :
     * `boltlmm` : 1 perform boltlmm (default 0), see boltlmm option in _plink-assoc.nf_
     * `gemma` : 1 perform gemma (default 0)  see gemma option in _plink-assoc.nf_
     * `covariates` : covariates to include in model (if ph_normalise is 1)
  * Statistics option :
     * `ph_alpha_lim` : list of alpha used to computed significance (separated by comma)  
     * `ph_windows_size` : windows size around position used to simulate phenotype to define if was detected, in bp ex 1000bp in CM ex 0.1CM

## output 
different output is provided :
   * simul folder : contains position used to defined phenotype 
   * in boltlmm/gemma folder,  res_boltlmm/gemma.stat  contains summary stat for each alpha:
      * we defined `windows true` as the windows around snp used to build phenotype (size is defined previously)
      * `nsig_simall_alpha` : number significant snp in all windows true 
      * `nsig_sim_alpha` :   number windows true where at least one snps is significant
      * `nsig_simaround_alpha` : number significant windows true where one snp is significant and has been excluded snps used to build pheno
      * `nsig_nosim_alpha` : snp significant snp not in windows true
      * `nsnp` : snp number total  in dataset
      * `nsnpsima` : snp number used to build phenotype (see ph_nb_qtl)
   * in boltlmm/gemma/simul/ : contains p.value compute for each simulation

#Note 
  * for phenotype simulation all missing values is discarded and replaced by more frequent allele
  * phenosim use a lot of memory and time, subsample of snp/samples improve times / memory used




# 11. Copyright and general

### Authors

Scott Hazelhurst, Lerato E. Magosi, Shaun Aron, Rob Clucas, Eugene de Beste, Aboyomini Mosaku, Don Armstrong and the Wits Bioinformatics team

We thank Harry Noyes from the University of Liverpool and Ayton Meintjes from UCT who both spent significant effort being testers of the pipleine.

### License
h3agwas offered under the MIT license. See LICENSE.txt.

### Download

`git clone https://github.com/h3abionet/h3agwas`

### References
