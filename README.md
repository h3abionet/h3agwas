
<img src="auxfiles/H3ABioNetlogo2.jpg"/>

# H3Agwas Pipeline Version 3

The major change from Version 2 to Version 3 is the reorganisation of the repo so that the different workflows are in separate directories.


This means that instead of running `nextflow run h3abionet/h3agwas/assoc.nf`, you should run `nextflow run h3abionet/h3agwas/assoc/main.nf`


## Brief introduction

In addition to this README we have a detailed tutorial and videos 
* These can be found at http://www.bioinf.wits.ac.za/gwas

pipeline do different step of GWAS :
 * [Format input illuminat in plink format](call2plink/README.md)
 * [Quality control of array input illuminat in plink format](qc/README.md)
 * [Association using different software : gcta, plink, gemma, Bolt-LMM, FastLMM and GxE with gemma and plink](assoc/README.md)
 * Post meta analyses script :
   * [meta analyse script and mtag approach](meta/README.md)
   * [Computation of heritabilities or variance explained of phenotype](heritabilities/README.md)
   * [Finemapping and cojo extraction of windows ](finemapping/README.md)
 * [Simulation of dataset](utils/build_example_data/README.md)
 * [Format data differents dataset](formatdata//README.md) :
   * plink in vcf to prepared your data at imputation 
   * vcf in plink after imputation


## What's new :
 * see [What's news](News.md)

## Background


H3Agwas is a simple human GWAS analysis workflow for data quality control (QC) and basic association testing developed by [H3ABioNet](https://www.h3abionet.org/). It is an extension of the [witsGWAS pipeline](http://magosil86.github.io/witsGWAS/) for human genome-wide association studies built at the [Sydney Brenner Institute for Molecular Bioscience](https://www.wits.ac.za/research/sbimb/). H3Agwas uses Nextflow as the basis for workflow managment and has been dockerised to facilitate portability.


The original version of the H3Agwas was published in June 2017 with minor updates and bug fixes through the rest of the year. Based on experience with large data sets, the pipelines were considerably revised with additional features, reporting and a slightly different workflow.  


We have moved all scripts from Python 2 to Python 3, so you will need to have Python 3 installed.  

_Please ignore the Wiki in this version which refers to version 1_


## Questions and feeback

Problems with the workflow should be raised as an issue on this  GitHub repo. (If you think the probem is the workflow)

If you need help with using the workflow, please log a call with the [H3A Help Disk](https://www.h3abionet.org/categories/communications/helpdesk)


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

There are three separate workflows that make up *h3agwas*

1. `call2plink`.  Conversion of Illumina genotyping reports with TOP/BOTTOM or FORWARD/REVERSE  calls into PLINK format, aligning the calls.
   * see [README of call2plink](call2plink/)

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

3. `assoc`: Association study. A simple analysis association study is done. The purpose of this is to give users an introduction to their data. Real studies, particularly those of the H3A consortium will have to handle compex co-variates and particular population study. We encourage users of our pipeline to submit their analysis for the use of other scientists.
   * see [README of assoc/](assoc/)
  * Basic PLINK association tests, producing manhattan and qqplots
  * CMH association test - Association analysis, accounting for clusters
  * permutation testing
  * logistic regression
  * Efficient Mixed Model Association testing with gemma, boltlmm or fastlmm
  * Gene environment association with gemma or plink

4. `meta` : meta analyse or mtag :
    * `meta/meta-assoc.nf` : do meta analysis with summary statistics 
    * `meta/mtag.nf` : do mtag analysis with summary statistics 

5. `heritabilities`
    *  ̀heritabilities/esth2-assoc.nf` : estimate heritability and co-heritabilie with gcta, ldsc, gemma and bolt

6. `finemapping` :
    * `finemapping/main.nf` : performed meta analysis using different data set 
    * `finemapping/cojo-assoc.nf` : do Conditional & joint (COJO) analysis of GWAS summary statistics without individual-level genotype data with gcta

7. `utils/build_example_data` 
   * `utils/build_example_data/main.nf` : extract data set from vcf file and simulate dataset 
   * `utils/build_example_data/simul-assoc.nf` : simulation of phenotype using phenosim 

8. `utils/permutatation` 
  * `utils/permutation/permutation-assoc.nf`: do a permutation test to reevaluate p.value with gemma
9. `formatdata` : additional script to format data added some missing information etc...
  *  see [README of formatdata/](formatdata/)



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

### Nextflow

[Nextflow](https://www.nextflow.io/) is a workflow language designed at the Centre for Genomic Regulation, Barcelona. Although it is a general workflow language for science, it comes out of a bioinformmatics group and strongly supports bioinformatics. 

Our pipeline is built using Nextflow. However, users do not need to know anything about Nextflow. Obviously if you can do some programming you can customise and extend the pipelines, but you do not need to know Nextflow yourself. 

Nextlow is very easy to install and is highly portable. It supports partial execution and pipelines that scale.  Nextflow supports our worklow requirements very well.

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

If you want to run the `assoc` pipeline then you should install gemma,fastlmm if you are using those options.

## 2.5 Installing the workflow

There are two approaches: let Nextflow manage this for you; or download using Git. The former is easier; you need to use Git if you want to change the workflow

### 2.5.1 Managing using Nextflow

To download the workflow you can say

`nextflow pull h3abionet/h3agwas`

If we update the workflow, the next time you run it, you will get a warning message. You can do another pull to bring it up to date.

If you manage the workflow this way, you will run the scripts, as follows
* `nextflow run h3abionet/h3agwas/call2plink/main.nf ..... `
* `nextflow run h3abionet/h3agwas/qc/main.nf ..... `
* `nextflow run h3abionet/h3agwas/assoc/main.nf ..... `

### 2.5.2 Managing with Git

Change directory where you want to install the software and say

    git clone https://github.com/h3abionet/h3agwas.git

This will create a directory called _h3agwas_ with all the necesssary code.
If you manage the workflow this way, you will run the scripts this way:
* `nextflow run SOME-PATH/call2plink ..... `
* `nextflow run SOME-PATH/qc ..... `
* `nextflow run SOME-PATH/assoc ..... `

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

`nextflow run h3abionet/h3agwas/qc/main.nf --input_dir=s3://h3abionet/sample`

If you have downloaded the sample data and the directory with the sample data is a sub-directory of your working directory, you could just say: `nextflow run h3abionet/h3agwas/qc/main.nf --input_dir=s3://h3abionet/sample`


### 3.1.2 If you downloaded using Git

Change directory to the directory in which the workflow was downloaded

`nextflow run  qc`

## 3.2 Remarks

The workflow runs and output goes to the _output_ directory. In the
_sampleA.pdf_ file, a record of the analysis can be found.

In order, to run the workflow on another PLINK data set, say _mydata.{bed,bim,fam}_, say

`nextflow run  qc --input_pat mydata`

(or `nextflow run  h3abionet/h3agwas/qc --input_pat mydata` : **for simplicity for the rest of the tutorial we'll only present the one way of running the workflow -- you should use the method that is appropriate for you**)

If the data is another directory, and you want to the data to go elsehwere:

`nextflow run  qc --input_pat mydata --input_dir /data/project10/ --output_dir ~/results `

There are many other options that can be passed on the the command-line. Options can also be given in the _config_ file (explained below). We recommend putting options in the configuration file since these can be archived, which makes the workflow more portable

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

```nextflow run -c data1.config qc```


**This is highly recommended.** We recommend that you keep the `nextflow.config` file as static as possible, perhaps not even modifying it from the default config. Then  for any
 run or data set, have a much smaller config file that only specifies the changes you want made. The base `nextflow.config` file will typically contain config options that are best set by the H3Agwas developers (e.g., the names of the docker containers) or default GWAS options that are unlikely to change. In your separate config file, you will specify the run-specific options, such as data sets, directories or particular GWAS parameters you want. Both configuration files should be specified. For example, suppose I create a sub-directory within the directory where the nextflow file is (probably called h3agwas). Within the h3agwas directory I keep my nexflow.config file and the nextflow file itself. From the sub-directory, I run the workflow by saying:

```nextflow run  -c data1.config ../qc```

This will automatically use the `nextflow.config` file in either the current or parent directory. Note that the the config files are processed in order: if an option is set into two config files, the latter one takes precedence.


## 4.2 Creating an auxiliary nextflow .config file

There is a template of a nextflow.config file called aux.config.template. This is a read only file. Make a copy of it, call it _aux.config_ (or some suitable name).  This file contains all the options a user is likely to want to change. It does not specify options like the names of docker containers etc. Of course, you can if you wish modify the nextflow.config file, but we recommend against it. Your auxiliary file should supplement the nextflow.config file.

Then fill in the details in the config that are required for your run. These are expained in more detail below.


## 4.3 Specifying options

When you run the the scripts there are a number of different options that you might want to use. These options are specified by using  the `-flag` or `--flag` notation. The flags with a single hyphen (e.g. `-resume`) are standard Nextflow options applicable to all Nextflow scripts. The flags with a double hyphen (e.g., `--pi_hat`) are options that are specific to _our_ scripts.  *Take care not to mix this up as it's an easy error to make, and may cause silent errors to occur.*


Almost all the workflow options that are in the _nextflow.config_ file can also be passed on the command line and they will then override anything in the config like. For example

```nextflow run qc --cut_miss  0.04```

sets the maximim allowable per-SNP misisng to 4%. However, this should only be used when debugging and playing round. Rather, keep the options in the auxiliary config file that you save. By putting options on the command line you reduce reproducibility. (Using the parameters that change the mode of the running -- e.g. whether using docker or whether to produce a time line only affects time taken and auxiliary data rather than the substantive results).


## 4.4 Partial execution and resuming execution

Often a workflow may fail in the middle of execution because there's a problem with data (perhaps a typo in the name of a file), or you may want to run the workflow with slightly different parameters. Nextflow is very good in detecting what parts of the workflow need to re-executed -- use the `-resume` option. 

## 4.5 Cleaning up 

If you want to clean up your work directory, say `nextflow clean`.

## 4.6 Workflow overview, and timing

Nextflow provides [several options](https://www.nextflow.io/docs/latest/tracing.html) for visualising and tracing workflow. See the Nextflow documentation for details. Two of the options are:

* A nice graphic of a run of your workflow

    `nextflow run qc -with-dag quality-d.pdf`

* A timeline of your workflow and individual processes (produced as an html file).

    `nextflow run <pipeline name> -with-timeline time.html`

    This is useful for seeing how long different parts of your process took. Also useful is peak virtual memory used, which you may need to know if running on very large data to ensure you have a big enough machine and specify the right parmeters.




# 5. Running the workflow in different environments

In the  quick start we gave an overview of running our workflows in different environments. Here we go through all the options, in a little more detail


## 5.1 Running natively on a machine

This option requires that all dependancies have been installed. You run the code by saying

````
nextflow run qc
````

You can add that any extra parameters at the end.


## 5.2 Running  on a local machine with Docker

This requires the user to have docker installed.

Run by `nextlow run qc -profile docker`


## 5.3 Running on a cluster 

Nextflow supports execution on clusters using standard resource managers, including Torque/PBS, SLURM and SGE. Log on to the head node of the cluster, and execute the workflow as shown below. Nextflow submits the jobs to the cluster on your behalf, taking care of any dependancies. If your job is likely to run for a long time because you've got really large data sets, use a tool like _screen_ to allow you to control your session without timing out.

Our workflow has pre-built configuration for SLURM and Torque/PBS. If you use another scheduler that Nextflow supports you'll need to do a _little_ more (see later): see https://www.nextflow.io/docs/latest/executor.html for details

To run using Torque/PBS, log into the head node. Edit the _nextflow.config_ file, and change the `queue` variable to be the queue you will run jobs on (if you're not sure of this, ask your friendly sysadmin). Then when you run, our workflow, use the `-profile pbs` option -- typically you would say something like `nextflow run -c my.config qc -profile pbs`. Note that the `-profile pbs` only uses a single "-".

Similarily, if you run SLURM, set the _queue_ variable, and use the `-profile slurm` option.

To use only of the other schedulers supported by Nextflow, add the following sub-stanza to your nextflow.config file inside of the _profile_ stanza:



```
    myscheduler {
        process.executor = 'myscheduler'
	process.queue = queue
    }
```

where `myscheduler` is one of: nqsii, htcondor, sge, lsf.  


and  then use this as the profile.


We assume all the data is visible to all nodes in the swarm. Log into the head node of the Swarm and run your chosed workflow -- for example


## 5.4 Running on Docker Swarm

We have tested our workflow on different Docker Swarms. How to set up Docker Swarm is beyond the scope of this tutorial, but if you have a Docker Swarm, it is easy to run. From the head node of your Docker swarm, run

```
nextflow run qc -profile dockerSwarm
```

## 5.5 Singularity


Our workflows now run easily with Singularity.

`nextflow run qc -profile singularity`

or

`nextflow run qc -profile pbsSingularity`

By default the user's ${HOME}/.singularity will be used as the cache for Singularity images. If you want to use something else, change the `singularity.cacheDir` parameter in the config file.


If you have a cluster which runs Docker, you can get the best of both worlds by editing the queue variable in the _pbsDocker_ stanza, and then running

```
nextflow run qc -profile option
```

where _option_ is one of _pbsDocker_, _pbsSingularity_, _slurmDocker_ or _slurmSingularity_. If you use a different scheduler, read the Nextflow documentation on schedulers, and then use what we have in the _nextflow.config_ file as a template to tweak.

## 5.5 Other container services

We are unlikely to support udocker unless Nextflow does. See this link for a discussion https://www.nextflow.io/blog/2016/more-fun-containers-hpc.html

## 5.6 Running on Amazon EC2


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


   __ Need to change : By default, running the workflow like this runs the `qc` script. If you want to run one of the other scripts you would say `nextflow run  h3abionet/h3agwas/topbottom.nf` or `nextflow run h3abionet/h3agwas/assoc.nf` etc. __


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

Note there are two uses of `-c`. The positions of these arguments are crucial. The first are arguments to _nextflow_ itself and gives the configuration files that nextflow to use. The second is an argument to _cloud create_ which says how many nodes should be created. 

The _scott.aws_ file is not shared or put under git control. The _nextflow.config_ and _run10.config_ files can be archived, put under git control and so on because you _want_ to share and archive this information with o thers.

## 5.7 Running on AWS Batch

AWS Batch is a service layered on top of EC2 by Amazon which may make it easier and / or cheaper than using EC2. My personal view is that if you are only our pipeline on Amazon and you have reasonable Linux experience then the EC2 implementation above is probably easier. However, if you use or plan to use AWS Batch for other services then, AWS Batch is a definite option.


### Step 1

Create an AWS Batch queue and computing environment.  Setting up AWS
Batch is beyond the scope of this document. You can look at [Amazon's
documentation](https://docs.aws.amazon.com/batch/latest/userguide/create-job-queue.html)
or the general documentation from BioNet.

You also need to set up an S3 bucket for working space. Remember to set permissions on this bucket appropriately.

### Step 2

Create a nextflow config file with your personal information (this should not be put under git !). Set the `process.queue` to the name of the queue you created in the previous step and replace the `accessKey`, `secretKey` and `region` parameters with your values.

```

process.queue = 'queue_name'

aws {
    accessKey ='accessKey'
    secretKey = 'WHATEVERYOURSECRETKEYISGOESHERE'
    region    ='eu-west-1'
}

```

You can call your config file whatever you want, but for sake of the documentation below I'm assuming you called in `aws.config`.

### Step 3



Set up your other config files as required. Note that data you wish to process  can either be local or in an S3 bucket.

### Step 4

Run the job (in this example the _qc_ worfklow).  You need to specify the s3 bucket to be used and also the `awsbatch` profile

```
nextflow run -c aws.config -c job.config qc  -bucket-dir s3://my-bucket/some/path  -profile awsbatch
```

# 6. Dealing with errors

One problem with our current workflow is that error messages can be obscure. Errors can be caused by
* bugs in our code
* your doing something odd

There are two related problems. When a Nextflow script fails for some reason, Nextflow prints out in _great_ detail what went wrong. Second, we don't always catch mistakes that the user makes gracefully.

First, don't panic. Take a breath and read through the error message to see if you can find a sensible error message there. 

A typical error message looks something like this

```
Command exit status:
  1

Command output:
  (empty)

Command error:
  Traceback (most recent call last):
    File ".command.sh", line 577, in <module>
      bfrm, btext = getBatchAnalysis()
    File ".command.sh", line 550, in getBatchAnalysis
      result = miss_vals(ifrm,bfrm,args.batch_col,args.sexcheck_report)
    File ".command.sh", line 188, in miss_vals
      g  = pd.merge(pfrm,ifrm,left_index=True,right_index=True,how='inner').groupby(pheno_col)
    File "/usr/local/python36/lib/python3.6/site-packages/pandas/core/generic.py", line 5162, in groupby
      **kwargs)
    File "/usr/local/python36/lib/python3.6/site-packages/pandas/core/groupby.py", line 1848, in groupby
      return klass(obj, by, **kwds)
    File "/usr/local/python36/lib/python3.6/site-packages/pandas/core/groupby.py", line 516, in __init__
      mutated=self.mutated)
    File "/usr/local/python36/lib/python3.6/site-packages/pandas/core/groupby.py", line 2934, in _get_grouper
      raise KeyError(gpr)

Column 'batches' unknown

Work dir:
  /project/h3abionet/h3agwas/test/work/cf/335b6d21ad75841e1e806178933d3d

Tip: when you have fixed the problem you can continue the execution appending to the nextflow command line the option `-resume`

 -- Check '.nextflow.log' file for details
WARN: Killing pending tasks (1)

```

Buried in this is an error message that might help (did you say there was a column _batches_ in the manifest?) If you're comfortable, you can change directory to the specified directory and explore. There'll you find
* Any input files for the process that failed
* Any output files that might have been created
* The script that was executed can be found in `.command.sh`
* Output and error can be found as `.command.out` and `.command.err`

If you spot the error, you can re-run the workflow (from the original directory), appending `-resume`.  Nextflow will re-run your workflow as needed -- any steps that finished successfully will not need to be re-run.

If you are still stuck you can ask for help at two places


* H3ABioNet Help desk --- https://www.h3abionet.org/support


* On GitHub -- need a GitHub account if you have a GitHub account

   https://github.com/h3abionet/h3agwas/issues


# 7. Auxiliary Programs

These are in the aux directory

## 7.1 updateFam.py

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

## 7.2 getRunsTimes.pl (By Harry Noyes)

Nextflow has great options for showing resourc usage. However, you have to remember to set those option when you run.  It's easy to forget to do this. This very useful script by Harry Noyes (harry@liverpool.ac.uk) parses the .nextflow.log file  for you

## 7.3 make_ref.py 

Makes a reference genome in a format the the pipeline can use. The first argument is a directory that contains FASTA files for each chromosome; the second is the strand report, the third is the manifest report, the fourt in the base of othe output files.


`python3 make_ref.py auxfiles/37/ H3Africa_2017_20021485_A3_StrandReport_FT.txt H3Africa_2017_20021485_A3.csv h3aref201812`


The program checks each SNP given in the manifest file by the chromosome and position number and then checks that the probe given in the manifest file actually matches the reference genome at that point. Minor slippage is acceptable because of indels.

The wrn file are SNPs which are probably OK but have high slippage (these are in the ref file)
The err file are the SNPs which don't match.

## 7.4 plates.py

This is used to depict where on the plates particular samples are. This is very useful for looking at problems in the data. If for example you find a bunch of sex mismatches this is most likely due to misplating. This script is a quick way of looking at the problem and seeing whether the errors are close together or spread out. There are two input arguments

* A file with the IDs of the individuals -- assuming that the first token on each line is an individual
* A sample sheet that gives the plating of each sample

There is one output parameter -- the name of a directory where output should go. The directory should exist.

You may need to change this line

```
batches['ID'] = batches['Institute Sample Label'].apply(lambda x:x[18:])
```

In our example, we assumed the ID can found in the column "Institute Sample Label" but from the position 18 (indexed from 0) in the string. Change as appropriate for you


# 8. Acknowledgement, Copyright and general

## Citing this workflow

If you use this workflow, please cite the following paper

* Baichoo S, Souilmi Y, Panji S, Botha G, Meintjes A, Hazelhurst S, Bendou H, De Beste E, Mpangase P, Souiai O, Alghali M, Yi L, O'Connor B, Crusoe M, Armstrong D, Aron S, Joubert D, Ahmed A, Mbiyavanga M, Van Heusden P, Magosi, L, Zermeno, J, Mainzer L, Fadlelmola F, Jongeneel CV, and Mulder N. (2018) Developing reproducible bioinformatics analysis workflows for heterogenous computing environments to support African genomics, *BMC Bioinformatics* **19**, 457, 13 pages, doi:10.1186/s12859-018-2446-1.


## Acknowledgement

We acknowledge funding by the National Institutes of Health through the NHGRI (U41HG006941). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

* We thank Sumir Panji and Nicola Mulder for their support and leadership
* We thank Fourie Joubert at the University of Pretoria for hosting our initial hackathon.


### Authors

Current team: Scott Hazelhurst, Jean-Tristan Brandenburg, Lindsay Clark, Obokula Smile, Michael Ebo Turkson, Michael Thompson, 

H3ABioNet Pipelines team leadership: Christopher Fields,  Shakuntala Baichoo, Sumir Panji, Gerrit Botha.


Past members and contributors: Lerato E. Magosi, Shaun Aron, Rob Clucas,  Eugene de Beste, Aboyomini Mosaku, Don Armstrong and the Wits Bioinformatics team

We thank Harry Noyes from the University of Liverpool and Ayton Meintjes from UCT who both spent significant effort being testers of the pipleine, and the many users at the Sydney Brenner Institute for Molecular Bioscience for their patience and suggestion.

### Licence
This software is licensed under the MIT Licence.

### Funding

We acknowledge the support from the NIH NHGRI H3ABioNet (U24HG006941)   and AWI-Gen   (U54HG006938)

### Download

`git clone https://github.com/h3abionet/h3agwas`


