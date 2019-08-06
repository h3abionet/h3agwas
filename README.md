


<img src="aux/H3ABioNetlogo2.jpg"/>

# h3agwas Pipeline Version 3

The major change from Version 2 to Version 3 is the reorganisation of the repo so that the different workflows are in separate directories.

This means that instead of running `nextflow run h3abionet/h3agwas/assoc.nf`, you should run `nextflow run h3abionet/h3agwas/assoc`


## Background



h3aGWAS is a simple human GWAS analysis workflow for data quality control (QC) and basic association testing developed by [H3ABioNet](https://www.h3abionet.org/). It is an extension of the [witsGWAS pipeline](http://magosil86.github.io/witsGWAS/) for human genome-wide association studies built at the [Sydney Brenner Institute for Molecular Bioscience](https://www.wits.ac.za/research/sbimb/). h3aGWAS uses Nextflow as the basis for workflow managment and has been dockerised to facilitate portability.



The original version of the h3aGWAS was published in June 2017 with minor updates and bug fixes through the rest of the year. Based on experience with large data sets, the pipelines were considerably revised with additional features, reporting and a slightly different workflow.  


We have moved all scripts from Python 2 to Python 3, so you will need to have Python 3 installed.  

_Please ignore the Wiki in this version which refers to version 1_


## Brief introduction

In addition to this README we have the following material available
* A short video overview of the pipeline can be found at http://www.bioinf.wits.ac.za/gwas/h3agwas.mp4
* A handout from a lecture can be found at http://www.bioinf.wits.ac.za/gwas/gwas-comp-handout.pdf

### Restrictions
This version has been run on real data sets and works. However, not all cases have been thoroughly tested. In particular
* it is not robust when X chromosome data is not available
* the reporting assumes you want to do batch/site analysis. If you don't the code works but the report may look a bit odd with some figures repeated.



## Outline of documentation

1. Features
2. Installing the pipeline
3. A quick start example
4. The Nextflow configuration file
5. The QC pipeline: `qc.nf`
6. A simple association testing pipeline: `assoc.nf`
7. Converting Illumina genotyping reports to PLINK: `topbottom.nf`
8. Advanced options: Docker, PBS, Singularity, Amazon EC2
9. Dealing with errors
10. Auxiliary Programs

#  Features

##  Goals of the h3agwas pipeline

The goals of this pipeline is to have a portable and robust pipeline
for performing a genome-wide association study

There are three separate workflows that make up *h3agwas*

1. `call2plink`.  Conversion of Illumina genotyping reports with TOP/BOTTOM or FORWARD/REVERSE  calls into PLINK format, aligning the calls.

2. `qc`: Quality control of the data. This is the focus of the pipeline. It takes as input PLINK data and has the following functions

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
  * Basic PLINK association tests, producing manhattan and qqplots
  * CMH association test - Association analysis, accounting for clusters
  * permutation testing
  * logistic regression
  * Efficient Mixed Model Association testing with gemma, boltlmm or fastlmm
  * Gene environment association with gemma or plink
  * Other scripts gave for post analysis :
    * `assoc/cojo-assoc.nf` : do Conditional & joint (COJO) analysis of GWAS summary statistics without individual-level genotype data with gcta
    * ̀ assoc/esth2-assoc.nf` : estimate heritability and co-heritabilie with gcta, ldsc, gemma and bolt
    * `assoc/meta-assoc.nf` : do meta analysis with summary statistics 
    * `assoc/permutation-assoc.nf`: do a permutation test to reevaluate p.value with gemma
    * `assoc/simul-assoc.nf` : simulation of bed file 




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

# 2. Installing h3aGWAS

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
* Java 8
* Nextflow. To install Nextflow, run the command below. It creates a _nextflow_ executable in the directory you ran the command. Move the executable to a directory on the system or user PATH and make it executable. You need to be running Nextflow 27 (January 2018) or later.
    `curl -fsSL get.nextflow.io | bash`

  If you don't have curl (you can use wget)

* Git 

## 2.3 Installing with Docker or Singularity

If you install Docker or Singularity, you do not need to install all the other dependencies. Docker is available on most major platforms.  See [the Docker documentation](https://docs.docker.com/) for installation for your platform.  Singularity works very well on Linux.
    
That's it. 

## 2.4  Installing software dependencies to run natively

This requires a standard Linux installation or macOS. It requires _bash_ to be available as the shell of the user running the pipeline.

The following code needs to be installed and placed in a directory on the user's PATH.

* plink 1.9 [Currently, it will not work on plink 2, though it is on our list of things to fix. It probably will work on plink 1.05 but just use plink 1.0]
* LaTeX. A standard installation of texlive should have all the packages you need. If you are installing a lightweight TeX version, you need the following pacakges which are part of texlive.: fancyhdr, datetime, geometry, graphicx, subfig, listings, longtable, array, booktabs, float, url.
* python 3.6 or later. pandas, numpy, scipy, matplotlib and openpyxl need to be installed. You can instally these by saying: `pip3 install pandas`  etc

If you want to run the `assoc.nf` pipeline then you should install gemma if you are using those options.

## 2.5 Installing the workflow

There are two approaches: let Nextflow manage this for you; or download using Git. The former is easier; you need to use Git if you want to change the workflow

### 2.5.1 Managing using Nextflow

To download the workflow you can say

`nextflow pull h3abionet/h3agwas`

If we update the workflow, the next time you run it, you will get a warning message. You can do another pull to bring it up to date.

If you manage the workflow this way, you will run the scripts, as follows
* `nextflow run h3abionet/h3agwas/topbottom.nf ..... `
* `nextflow run h3abionet/h3agwas/qc.nf ..... `
* `nextflow run h3abionet/h3agwas/assoc.nf ..... `

### 2.5.2 Managing with Git

Change directory where you want to install the software and say

    `git clone https://github.com/h3agwas`

This will create a directory called _h3agwas_ with all the necesssary code.
If you manage the workflow this way, you will run the scripts this way:
* `nextflow run SOME-PATH/topbottom.nf ..... `
* `nextflow run SOME-PATH/qc.nf ..... `
* `nextflow run SOME-PATH/assoc.nf ..... `

where _SOME-PATH_ is a relative or absolute path to where the workflow was downloaded.



# 3. Quick start example

This section shows a simple run of the `qc` pipeline that
should run out of the box if you have installed the software or
Docker. More details and general configuration will be shown later.

This section illustrates how to run the pipeline on a small sample data
file with default parameters.  For real runs, the data to be analysed
and the various parameters to be used are specified in the
_nextflow.config_ file.  The details will be explained in another
section.

If you have downloaded the software using Git, you can find the sample data in the directory. Otherwise you can download the files from http://www.bioinf.wits.ac.za/gwas/sample.sip and unzip


The sample data to be used is in the _input_ directory (in PLINK
format as _sampleA.bed_, _sampleA.bim_, _sampleA.fam_). The default
_nextflow.config_ file uses this, and so you can run the workflow
through with this example. Note that this is a very small PLINK data set 
with no X-chromosome information and no sex checking is done.




## 3.1 Running on your local computer 

This requires that all software dependancies have been installed. 

### 3.1.1 If you downloaded using Nextflow

We also assume the _sample_ directory with data is in the current working directory

`nextflow run h3abionet/h3agwas/qc`


### 3.1.2 If you downloaded using Git

Change directory to the directory in which the workflow was downloaded

`nextflow run  qc`

## 3.2 Remarks

The workflow runs and output goes to the _output_ directory. In the
_sampleA.pdf_ file, a record of the analysis can be found.

In order, to run the workflow on another PLINK data set, say _mydata.{bed,bim,fam}_, say

`nextflow run  qc --input_pat mydata`

(or `nextflow run  h3abionet/h3agwas/qc.nf --input_pat mydata` : **for simplicity for the rest of the tutorial we'll only present the one way of running the workflow -- you should use the method that is appropriate for you**)

If the data is another directory, and you want to the data to go elsehwere:

`nextflow run  qc.nf --input_pat mydata --input_dir /data/project10/ --output_dir ~/results `

There are many other options that can be passed on the the command-line. Options can also be given in the _config_ file (explained below). We recommend putting options in the configuration file since these can be archived, which makes the workflow more portable

## 3.3 Running with Docker on your local computer

Execute 

`nextflow run  qc.nf -profile docker`

Please note that the _first_ time you run the workflow using Docker,  the Docker images will be downloaded. *Warning:* This will take about 1GB of bandwidth which will consume bandwidth and will take time depending on your network connection. It is only the first time that the workflow runs that the image will be downloaded.


More options are shown later.

##3.4 Running multiple workflows at the same time

You may at some point want to run multiple, _independent_ workflows at the same time (e.g. different data). This is possible. However, each run should be started in a different working directory. You can refer to the scripts and even the data in the same diretory, but the directories from which you run the `nextflow run` command should be different.


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

```nextflow run -c data1.config qc.nf```


**This is highly recommended.** We recommend that you keep the `nextflow.config` file as static as possible, perhaps not even modifying it from the default config. Then  for any
 run or data set, have a much smaller config file that only specifies the changes you want made. The base `nextflow.config` file will typically contain config options that are best set by the h3aGWAS developers (e.g., the names of the docker containers) or default GWAS options that are unlikely to change. In your separate config file, you will specify the run-specific options, such as data sets, directories or particular GWAS parameters you want. Both configuration files should be specified. For example, suppose I create a sub-directory within the directory where the nextflow file is (probably called h3agwas). Within the h3agwas directory I keep my nexflow.config file and the nextflow file itself. From the sub-directory, I run the workflow by saying:

```nextflow run  -c data1.config ../qc.nf```

This will automatically use the `nextflow.config` file in either the current or parent directory. Note that the the config files are processed in order: if an option is set into two config files, the latter one takes precedence.


## 4.2 Creating an auxiliary nextflow .config file

There is a template of a nextflow.config file called aux.config.template. This is a read only file. Make a copy of it, call it _aux.config_ (or some suitable name).  This file contains all the options a user is likely to want to change. It does not specify options like the names of docker containers etc. Of course, you can if you wish modify the nextflow.config file, but we recommend against it. Your auxiliary file should supplement the nextflow.config file.

Then fill in the details in the config that are required for your run. These are expained in more detail below.


## 4.3 Specifying options

When you run the the scripts there are a number of different options that you might want to use. These options are specified by using  the `-flag` or `--flag` notation. The flags with a single hyphen (e.g. `-resume`) are standard Nextflow options applicable to all Nextflow scripts. The flags with a double hyphen (e.g., `--pi_hat`) are options that are specific to _our_ scripts.  *Take care not to mix this up as it's an easy error to make, and may cause silent errors to occur.*


Almost all the workflow options that are in the _nextflow.config_ file can also be passed on the command line and they will then override anything in the config like. For example

```nextflow run qc.nf   --cut_miss  0.04```

sets the maximim allowable per-SNP misisng to 4%. However, this should only be used when debugging and playing round. Rather, keep the options in the auxiliary config file that you save. By putting options on the command line you reduce reproducibility. (Using the parameters that change the mode of the running -- e.g. whether using docker or whether to produce a time line only affects time taken and auxiliary data rather than the substantive results).


## 4.4 Partial execution and resuming execution

Often a workflow may fail in the middle of execution because there's a problem with data (perhaps a typo in the name of a file), or you may want to run the workflow with slightly different parameters. Nextflow is very good in detecting what parts of the workflow need to re-executed -- use the `-resume` option. 

## 4.5 Cleaning up 

If you want to clean up your work directory, say `nextflow clean`.

## 4.6 Workflow overview, and timing

Nextflow provides [several options](https://www.nextflow.io/docs/latest/tracing.html) for visualising and tracing workflow. See the Nextflow documentation for details. Two of the options are:

* A nice graphic of a run of your workflow

    `nextflow run qc.nf -with-dag quality-d.pdf`

* A timeline of your workflow and individual processes (produced as an html file).

    `nextflow run <pipeline name> -with-timeline time.html`

    This is useful for seeing how long different parts of your process took. Also useful is peak virtual memory used, which you may need to know if running on very large data to ensure you have a big enough machine and specify the right parmeters.









# 5. Running the workflow in different environments

In the  quick start we gave an overview of running our workflows in different environments. Here we go through all the options, in a little more detail

## 5.1 Running natively on a machine

This option requires that all dependancies have been installed. You run the code by saying

````
nextflow run qc.nf
````

You can add that any extra parameters at the end.

## 5.2 Running  on a local machine with Docker

This requires the user to have docker installed.

Run by `nextlow run qc.nf -profile docker`



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
nextflow run qc.nf -profile dockerSwarm
```

## 5.5 Singularity


Our workflows now run easily with Singularity.

`nextflow run qc.nf -profile singularity`

or

`nextflow run qc.nf -profile pbsSingularity`

By default the user's ${HOME}/.singularity will be used as the cache for Singularity images. If you want to use something else, change the `singularity.cacheDir` parameter in the config file.

## Running on a cluster with Docker or Singularity

If you have a cluster which runs Docker, you can get the best of both worlds by editing the queue variable in the _pbsDocker_ stanza, and then running

```
nextflow run qc.nf -profile option
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


   By default, running the workflow like this runs the `qc.nf` script. If you want to run one of the other scripts you would say `nextflow run  h3abionet/h3agwas/topbottom.nf` or `nextflow run h3abionet/h3agwas/assoc.nf` etc. 


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


#9. Dealing with errors

One problem with our current workflow is that error messages can be obscure. Errors can be caused by
* bugs in our code
* you doing something odd

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


# 6. Auxiliary Programs

These are in the aux directory

## 6.1 updateFam.py

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

## 6.2 getRunsTimes.pl (By Harry Noyes)

Nextflow has great options for showing resourc usage. However, you have to remember to set those option when you run.  It's easy to forget to do this. This very useful script by Harry Noyes (harry@liverpool.ac.uk) parses the .nextflow.log file  for you

## 6.3 make_ref.py 

Makes a reference genome in a format the the pipeline can use. The first argument is a directory that contains FASTA files for each chromosome; the second is the strand report, the third is the manifest report, the fourt in the base of othe output files.


`python3 make_ref.py aux/37/ H3Africa_2017_20021485_A3_StrandReport_FT.txt H3Africa_2017_20021485_A3.csv h3aref201812`


The program checks each SNP given in the manifest file by the chromosome and position number and then checks that the probe given in the manifest file actually matches the reference genome at that point. Minor slippage is acceptable because of indels.

The wrn file are SNPs which are probably OK but have high slippage (these are in the ref file)
The err file are the SNPs which don't match.

## 6.6 plates.py

This is used to depict where on the plates particular samples are. This is very useful for looking at problems in the data. If for example you find a bunch of sex mismatches this is most likely due to misplating. This script is a quick way of looking at the problem and seeing whether the errors are close together or spread out. There are two input arguments

* A file with the IDs of the individuals -- assuming that the first token on each line is an individual
* A sample sheet that gives the plating of each sample

There is one output parameter -- the name of a directory where output should go. The directory should exist.

You may need to change this line

```
batches['ID'] = batches['Institute Sample Label'].apply(lambda x:x[18:])
```

In our example, we assumed the ID can found in the column "Institute Sample Label" but from the position 18 (indexed from 0) in the string. Change as appropriate for you



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
  * `num_cores` : cores number used 
  * `ph_mem_req` : memory request for phenosim
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
     * `boltlmm` : 1 perform boltlmm (default 0), see boltlmm option in _assoc.nf_
     * `gemma` : 1 perform gemma (default 0)  see gemma option in _assoc.nf_
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




# 11 MetaAnalysis pipeline : `meta-assoc.nf`

This section describes a pipeline in devlopment, purpose of this pipeline is to do a meta analysis with a various format files.Our script, *meta-assoc.nf* takes as input various GWAS results files and `rsid` to do a metanalysis with METAL, GWAMA and Metasoft

## Installation
need python3, METAL, GWAMA, MR-MEGA and MetaSoft

## Running
The pipeline is run: `nextflow run meta-assoc.nf`


The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * `input`, `output` and script directories: the default is that these are subdirectories of the `work_dir` and there'll seldom be reason to change these;
  * `output_dir` = "all"
  * meta analysis option :
     * `metal` : 1 perform metal (default 0) 
     * `gwama` : 1 perform gwama (default 0)
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
       * `metal_bin` : binarie for metal (default : _metal_ ) 
       * `gwama_bin` :  binarie for gwam ( default : _GWAMA__ )
       * `metasoft_bin` : binarie for java of metasoft ( default _Metasoft.jar_)
       * `mrmega_bin` : binarie for java of metasoft ( default _Metasoft.jar_)
     * options softwares :
       * `ma_metasoft_opt` : append other option in metasoft command line(default : null)
       * `ma_genomic_cont` : use a genomic_control use in METAL and GWAMA(default, 0)
       * `ma_inv_var_weigth`: do a invert variance weight usefull for metal (default, 0)
       * `ma_random_effect` : do mixed model (default 1)
       * `ma_mrmega_pc` : how many pcs used for mrmega (default : 4)
       * `ma_mrmega_opt` : append other option in MR-MEGA command line (default : null)
## specificity 
### MR-MEGA
MR-MEGA need chromosomes, positions and N (sample number) for each position, so in pipeline referent file (in file_config, 1 in IsRefFile) must be have chromosome and poosition 


# 12 annotation pipeline: `annot-assoc.nf`
This section describes a pipeline in devlopment, objectives is annotation of rs using annotation, locuszoom, and phenotype in function of genotype

## Installation
need locuszoom, _R_ : (ggplot2), python3

## Running
The pipeline is run: `nextflow run annot-assoc.nf`

The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * input, output and script directories: the default is that these are subdirectories of the `work_dir` and there'll seldom be reason to change these;
  * `input_pat` : this typically will be the base name of the PLINK files you want to process (i.e., do not include the file suffix). But you could be put any Unix-style glob here. The workflow will match files in the relevant `input_dir` directory;

# 13. COJO : `cojo-assoc.nf`
this section describes a pipeline in devloment, objectives is doing a conditional and joint association using GWAS summary data and gcta
see [cojo](https://cnsgenomics.com/software/gcta/#COJO)
## Installation
need python3, gcta
## Running
The pipeline is run: `nextflow run annot-assoc.nf`

The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * `output_dir` : output directory
  * `output` : output pattern 
  * ̀ data` : same option that _plink-assoc.nf_, file is optional, used if need select specific individus for gcta,  compute frequencies or N, if mission in `file_gwas`
  * `input_pat`: the base of set of PLINK bed,bim and fam files (this should only match one);
  * `pheno` : optional, header in data, if present select individuals with no missiong individual to keep individuals for computed frequencie or gcta
  * `cut_maf` minor allele frequencies [ default : 0.0001]
  * ̀`file_gwas` : file contains gwas result, if N or frequencies is not available, it is computed with plink file and `data` file, to change format header must be defined :
    * ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
    * `head_freq` : freq header [ default : None], if not present computed with plink, (and data/pheno if present)
    * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
    * `head_rs` : rs header column [default : "SNP"]
    * `head_beta` : beta header colum [default : "BETA"]
    * `head_se`  : column for standard error of beta "SE"
    * `head_A1` : column for A0 :[default : "ALLELE0" ]
    * `head_A2` : column for A0 :[default : "ALLELE2" ]

Cojo parameter :
  * `cojo_wind` :  Specify a distance d (in Kb unit). It is assumed that SNPs more than d Kb away from each other are in complete linkage equilibrium. The default value is 10000 Kb (i.e. 10 Mb) if not specified. [ default : 10000 ]
  * `cojo_actual_geno` : If the individual-level genotype data of the discovery set are available (e.g. a single-cohort GWAS), you can use the discovery set as the reference sample. *option to avoid due to a various bug*  [default 0]
  * `cojo_slct` : Perform a stepwise model selection procedure to select independently associated SNPs? 1 : yes 0 : no [default 1] 
    * `cojo_p` :  Threshold p-value to declare a genome-wide significant hit. The default value is 5e-8 if not specified. This option is only valid in conjunction with the option `cojo_slct`. 
    * `cojo_slct_other` : other option for slct see [manual](https://cnsgenomics.com/software/gcta/#COJO)
  * `cojo_top_snps_chro` :  Perform a stepwise model selection procedure to select a fixed number of independently associated SNPs by chromosome without a p-value threshold.  [integer between 0 and n, to define top snp number. default : 0].
  * `gcta_mem_req`="6GB"
  
# 13. h2 estimation : `esth2-assoc.nf`

This section describes a pipeline in devlopment, objectives is estimated heritabilities with various way, we developped : ldlc, grmel of bolt and greml of gcta, gemma  
two distincs approaches should be considered :
  * based on relatdness matrix and phenotype as gcta, bolt, gemma  
  * based on gwas result as implemented in ldlc and gemma

## Installation
need python3, gcta, ldlc, bolt and gemma

## Running
The pipeline is run: `nextflow run esth2-assoc.nf`
 
The key options are:
  * `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
  * `output_dir` : output directory
  * `output` : output pattern
  * ̀ data` : same option that _plink-assoc.nf_, file is optional, used for gemma, bolt and gcta
    * `pheno` :phenotypes used in data to computed in gemma, bolt
  * `file_gwas` : one ore more one file gwas, used for ldsc and gemma, to defined column of files :
    * ̀ head_pval` : pvalue header [ default : "P_BOLT_LMM" ]
    * `head_n` : N (individuals number) [ default : None ], if not present computed with plink (and data/pheno if present)
    * `head_rs` : rs header column [default : "SNP"]
    * `head_beta` : beta header colum [default : "BETA"]
    * `head_se`  : column for standard error of beta "SE"
    * `head_A1` : column for A0 :[default : "ALLELE0" ]
    * `head_A2` : column for A0 :[default : "ALLELE2" ]
    * `head_freq` : freq header [ default : A1Freq], 
    * `head_n`: N header, used just for ldsc, if not present, `Nind` must be initialize.
    * `Nind` : if `head_n` not initialise, must be initialise, individuals number for each gwas file, separate by comma
  * `ldsc_h2` : need a estimation of h2 by ldc : 1 [default : 0]:   
    * [LDSC](https://github.com/bulik/ldsc) computes heritabilies between gwas value using LD information.
    * `ldsc_bin` : binary for ldsc 
    * `dir_ref_ld_chr` : folder containing ld information, 1 file by begin by chromosome num without chr, see : `--ref-ld-chr` in ldsc manual
    * ̀`ldsc_mem_req`
    * `ldsc_h2_multi` : computing genetic correlation. between different gwas result
    * `munge_sumstats_bin` : binary for munge 
    * `ldsc_h2opt` : other option for ldsc  
    * output : 
  * `gemma_h2` :
  * `gemma_h2_pval`
    * _-vc_ for vc equal 2 need LD files, fou




# 14. Acknowledgement, Copyright and general

## Acknowledgement

We acknowledge funding by the National Institutes of Health through the NHGRI (U41HG006941). The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

* We thank Sumir Panji and Nicola Mulder for their support and leadership
* We thank Fourie Joubert at the University of Pretoria for hosting our initial hackathon.
>>>>>>> master

### Authors

Scott Hazelhurst, Lerato E. Magosi, Shaun Aron, Rob Clucas, Jean-Tristan Brandenburg, Eugene de Beste, Aboyomini Mosaku, Don Armstrong and the Wits Bioinformatics team

We thank Harry Noyes from the University of Liverpool and Ayton Meintjes from UCT who both spent significant effort being testers of the pipleine.

### License
This software is licensed under the MIT Licence.


### Download

`git clone https://github.com/h3abionet/h3agwas`

### References
