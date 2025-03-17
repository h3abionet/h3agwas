<img src="helperfiles/H3ABioNetlogo2.jpg"/>

# H3Agwas Pipeline Version 4

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/5145/badge)](https://bestpractices.coreinfrastructure.org/projects/5145)
[![biotools:h3agwas](https://img.shields.io/badge/biotools-h3agwas-blue)](https://bio.tools/h3agwas)
[![Docker Cloud Build Status](https://img.shields.io/docker/cloud/build/h3abionet/h3agwas)](https://quay.io/organization/h3abionet_org)
[![Amazon aws](https://img.shields.io/badge/aws-pass-green)](Readme_AWS_Batch.md)
![version](https://img.shields.io/github/v/release/h3abionet/h3agwas)
[![DOI](https://zenodo.org/badge/66284716.svg)](https://zenodo.org/doi/10.5281/zenodo.3235520)

## How to cite

***Found publication in [BMC bio-informatics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05034-w/)***

***found last update in [dev branch](https://github.com/h3abionet/h3agwas/tree/dev)***

## Major change from previous version
 
The major change from Version 3 to Version 4 :
* dsl 2 supported 
* X, Y chromosomes supported for quality control, format dataset and will be integred
* imputation had been aded

## Brief introduction

In addition to this README we have a detailed tutorial and videos
* found example data-set and example in  [h3agwas-examples github](https://github.com/h3abionet/h3agwas-examples)

pipeline do different step of GWAS :
 * [Deleted technical duplicated](qc_dup/README.md)
 * [Quality control of array plink format](qc/README.md)
 *

## What's new :

 * Update see [What's new](News.md)
 * to see conversion in [V4](v4.md)

##  Outline of documentation
H3Agwas is a fully human GWAS analysis workflow for data quality control (QC), association and post association testing developed by [H3ABioNet](https://www.h3abionet.org/). It is an extension of the [witsGWAS pipeline](http://magosil86.github.io/witsGWAS/) for human genome-wide association studies built at the [Sydney Brenner Institute for Molecular Bioscience](https://www.wits.ac.za/research/sbimb/). H3Agwas uses Nextflow as the basis for workflow managment and has been dockerised to facilitate portability.


The original version of the H3Agwas was published in June 2017 with minor updates and bug fixes through the rest of the year. Based on experience with large data sets, the pipelines were considerably revised with additional features, reporting and a slightly different workflow.


We have moved all scripts from Python 2 to Python 3, so you will need to have Python 3 installed.

_Please ignore the Wiki in this version which refers to version 1_


## Questions and feedback

Problems with the workflow should be raised as an [issue](https://github.com/h3abionet/h3agwas/issues) on this  GitHub repo. (If you think the probem is the workflow) or contact [Dr. jean-tristan brandenburg](mailto:jean-tristan.brandenburg@wits.ac.za?subject=[H3AGWAS]%20questions)

If you need help with using the workflow, please log a call with the [H3A Help Desk](https://www.h3abionet.org/categories/communications/helpdesk) or contact [Dr. jean-tristan brandenburg](mailto:jean-tristan.brandenburg@wits.ac.za?subject=[H3AGWAS]%20questions)

If you need help with using the workflow, please log a call with the [H3A Help Desk](https://www.h3abionet.org/categories/communications/helpdesk) or contact [Dr. jean-tristan brandenburg](mailto:jean-tristan.brandenburg@wits.ac.za?subject=[H3AGWAS]%20questions)


## Outline of documentation
### main documentation
1. Features
2. Installing the pipeline
3. A quick start example
4. The Nextflow configuration file
5. Running the workflow in different environments and Advanced options: Docker, PBS, Singularity, Amazon EC2
6. Dealing with errors
7. Auxiliary Programs
8. Acknowledgement, Copyright and general

# 1.  Features

##  Goals of the h3agwas pipeline

The goals of this pipeline is to have a portable and robust pipeline
for performing a genome-wide association study

1. calling raw data from genotyping in vcf and plink
 - `` : 
 - `call2plink` Conversion of Illumina genotyping reports with TOP/BOTTOM or FORWARD/REVERSE  calls into PLINK format, aligning the calls.


2. `qc`: Quality control of the data. This is the focus of the pipeline. It takes as input PLINK data and has the following functions
   * see [README of qc](qc/)
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
3. `qc`
## Design principles

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


### Docker

A GWAS requires several software tools to be installed. Using Docker we can simplify the installation. Essentially, Docker wraps up all software dependancies into _containers_. Instead of installing all the dependancies, you can install Docker, easily and then install our containers. (In fact you don't need to explicitly install our containers, Nextflow and our workflow will do that for you automatically).

We expect that many of our users will use Docker. However, we recognise that this won't be suitable for everyone because many high performance computing centres do not support Docker for security reasons. It is possible to run our pipeline without Docker and will give  instructions about which software needs to be installed.

### Singularity

Similarily we support Singularity. Although it's a new feature, we've tested it two different organisaitons and it's worked flawlessly

# 2. Installing H3Agwas

## 2.1 Background

The h3agwas pipeline can be run in different environments; the requirements differ. The different modes are described in detail below
* Running on Docker/Singularity. This is the easiest way of running h3agwas. We have a set of Docker containers that have all the required executables and libraries.
* Running natively on a local computer -- this is requires a number of external executables and libraries to be installed..
* Running with a scheduler -- Nextflow supports a range of schedulers. Our pipeline supports using docker or running natively.
* Running on Amazon EC2.  You need to have Amazon AWS credentials (and a credit card). Our EC2 pipeline uses Docker so this is very easy to run.
* We have also used Docker swarm. If you have a Docker swarm it's easy to do.

We now explore these in details

## 2.2 Pre-requisites

**All** modes of h3agwas have the following requirements
* Java 8 or later
* Nextflow. To install Nextflow, run the command below. It creates a _nextflow_ executable in the directory you ran the command. Move the executable to a directory on the system or user PATH and make it executable. You need to be running Nextflow 27 (January 2018) or later.
    `curl -fsSL get.nextflow.io | bash`

  If you don't have curl (you can use wget)

* Git (this probably is already installed)

## 2.3 Installing with Docker or Singularity

If you install Docker or Singularity, you do not need to install all the other dependencies. Docker is available on most major platforms.  See [the Docker documentation](https://docs.docker.com/) for installation for your platform.  Singularity works very well on Linux.

That's it.

## 2.4  Installing software dependencies to run natively

This requires a standard Linux installation or macOS. It requires _bash_ to be available as the shell of the user running the pipeline.

The following code needs to be installed and placed in a directory on the user's PATH.

* plink 1.9 [Currently, it will not work on plink 2, though it is on our list of things to fix. It probably will work on plink 1.05 but just use plink 1.0]
* LaTeX. A standard installation of texlive should have all the packages you need. If you are installing a lightweight TeX version, you need the following pacakges which are part of texlive.: fancyhdr, datetime, geometry, graphicx, subfig, listings, longtable, array, booktabs, float, url.
* python 3.6 or later. pandas, numpy, scipy, matplotlib and openpyxl need to be installed. You can instally these by saying: `pip3 install pandas`  etc

If you want to run other scripts than qc then you should install various software of association, finemaping etc...

[you can find in h3agwas-example a dictatorial to installsoftwares on ubuntu platform](https://github.com/h3abionet/h3agwas-examples/tree/main/requirement/ubuntu)


## 2.5 Installing the workflow

There are two approaches: let Nextflow manage this for you; or download using Git. The former is easier; you need to use Git if you want to change the workflow

### 2.5.1 Managing using Nextflow

To download the workflow you can say

`nextflow pull h3abionet/h3agwas`

If we update the workflow, the next time you run it, you will get a warning message. You can do another pull to bring it up to date.

If you manage the workflow this way, you will run the scripts, as follows
* `nextflow run h3abionet/h3agwas/main.nf --qc 1  ..... `

### 2.5.2 Managing with Git

Change directory where you want to install the software and say

    git clone https://github.com/h3abionet/h3agwas.git

This will create a directory called _h3agwas_ with all the necesssary code.
If you manage the workflow this way, you will run the scripts this way:

where _SOME-PATH_ is a relative or absolute path to where the workflow was downloaded.

# 3. Quick start example

This section shows a simple run of the `qc` pipeline that
should run out of the box if you have installed the software or
Docker or Singularity. More details and general configuration will be shown later.

This section illustrates how to run the pipeline on a small sample data
file with default parameters.  For real runs, the data to be analysed
and the various parameters to be used are specified in the
_nextflow.config_ files in assoc, qc and call2plink folder.  The details will be explained in another
section.


Our quick start example will fetch the data from an Amazon S3 bucket,
but if you'd prefer then you use locally installed sample. If you have
downloaded the software using Git, you can find the sample data in the
directory. Otherwise you can download the files from
http://www.bioinf.wits.ac.za/gwas/sample.zip and unzip The sample data
to be used is in the _input_ directory (in PLINK format as
_sampleA.bed_, _sampleA.bim_, _sampleA.fam_). The default
_nextflow.config_ file uses this, and so you can run the workflow
through with this example. Note that this is a very small PLINK data
set with no X-chromosome information and no sex checking is done.


## 3.1 Running on your local computer

This requires that all software dependancies have been installed (see later for singularity or docker)

### 3.1.1 If you downloaded using Nextflow

We also assume the _sample_ directory with data is in the current working directory

`nextflow run h3abionet/h3agwas/main.nf --input_dir=s3://h3abionet/sample`

If you have downloaded the sample data and the directory with the sample data is a sub-directory of your working directory, you could just say: `nextflow run h3abionet/h3agwas/qc/main.nf --input_dir=s3://h3abionet/sample`

## 3.2 Remarks


## 3.3 Running with Docker on your local computer

Just add `-profile docker` to your run command -- for example,

* `nextflow run  qc -profile docker`   or
* `nextflow run h3abionet/h3agwas/qc/main.nf --input_dir=s3://h3abionet/sample`

Please note that the _first_ time you run the workflow using Docker,  the Docker images will be downloaded. *Warning:* This will take about 1GB of bandwidth which will consume bandwidth and will take time depending on your network connection. It is only the first time that the workflow runs that the image will be downloaded.


More options are shown later.

## 3.4 Running multiple workflows at the same time

You may at some point want to run multiple, _independent_ executions of the workflows at the same time (e.g. different data). This is possible. However, each run should be started in a different working directory. You can refer to the scripts and even the data in the same diretory, but the directories from which you run the `nextflow run` command should be different.


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

```nextflow run -c data1.config main.nf```


**This is highly recommended.** We recommend that you keep the `nextflow.config` file as static as possible, perhaps not even modifying it from the default config. Then  for any
 run or data set, have a much smaller config file that only specifies the changes you want made. The base `nextflow.config` file will typically contain config options that are best set by the H3Agwas developers (e.g., the names of the docker containers) or default GWAS options that are unlikely to change. In your separate config file, you will specify the run-specific options, such as data sets, directories or particular GWAS parameters you want. Both configuration files should be specified. For example, suppose I create a sub-directory within the directory where the nextflow file is (probably called h3agwas). Within the h3agwas directory I keep my nexflow.config file and the nextflow file itself. From the sub-directory, I run the workflow by saying:

```nextflow run  -c data1.config ../main.nf```

This will automatically use the `nextflow.config` file in either the current or parent directory. Note that the the config files are processed in order: if an option is set into two config files, the latter one takes precedence.

## 4.2 Creating an auxiliary nextflow .config file

There is a template of a nextflow.config file called aux.config.template. This is a read only file. Make a copy of it, call it _aux.config_ (or some suitable name).  This file contains all the options a user is likely to want to change. It does not specify options like the names of docker containers etc. Of course, you can if you wish modify the nextflow.config file, but we recommend against it. Your auxiliary file should supplement the nextflow.config file.

Then fill in the details in the config that are required for your run. These are expained in more detail below.

# 4.3 Specifying options

When you run the the scripts there are a number of different options that you might want to use. These options are specified by using  the `-flag` or `--flag` notation. The flags with a single hyphen (e.g. `-resume`) are standard Nextflow options applicable to all Nextflow scripts. The flags with a double hyphen (e.g., `--pi_hat`) are options that are specific to _our_ scripts.  *Take care not to mix this up as it's an easy error to make, and may cause silent errors to occur.*


Almost all the workflow options that are in the _nextflow.config_ file can also be passed on the command line and they will then override anything in the config like. For example

```nextflow run main.f --cut_miss  0.04```

sets the maximim allowable per-SNP misisng to 4%. However, this should only be used when debugging and playing round. Rather, keep the options in the auxiliary config file that you save. By putting options on the command line you reduce reproducibility. (Using the parameters that change the mode of the running -- e.g. whether using docker or whether to produce a time line only affects time taken and auxiliary data rather than the substantive results).

