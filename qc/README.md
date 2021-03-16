<img src="../auxfiles/H3ABioNetlogo2.jpg"/>

# The QC pipeline: `qc.nf`


This section describes the various ways in which the pipeline can be run and various options. Usually options are specified in the _nextflow.config_ file (or which ever file you use). However, you can also pass parameters to the Nextflow script on the command-line. Parameters on the command line over-ride any parameters specified in the config file. If you are a first time user of the pipeline you should read the discussion of the `nextflow.config` file in the parent README.


The main pipeline is the PLINK QC pipeline. It takes as input PLINK bed,bim,fam files and performs quality control on  the data according to the parameters specified in the config file.

The Nextflow script file is called *qc.nf*. This could be
called, for example, by running `nextflow run qc`.

The output of the QC is a set of PLINK files that can be used for GWAS, as well as PDF report that describes the QC steps.

### Restrictions
This version has been run on real data sets and works. However, not all cases have been thoroughly tested. In particular
* it is not robust when X chromosome data is not available
* the reporting assumes you want to do batch/site analysis. If you don't the code works but the report may look a bit odd with some figures repeated.

## 1 Input/Output :  PLINK format

Users will run the pipeline giving as input PLINK 1.9 bed, bim and fam files.  The key Nextflow parameters to set are:

* `work_dir` : the directory in which you will run the workflow. This will typically be the _h3agwas_ directory which you cloned;
* input, output and script directories: the default is that these are subdirectories of the `work_dir` and there'll seldom be reason to change these;
* `input_pat` : this typically will be the base name of the PLINK files you want to process (i.e., do not include the file suffix). But you could be put any Unix-style glob here. The workflow will match files in the relevant `input_dir` directory;
* `high_ld_regions_fname`: this is optional -- it is a list of regions which are in very high LD -- and are exclude when checking for relationships (https://www.cog-genomics.org/plink/1.9/filter#mrange_id).  Provide either absolute file path or relative to where you are running. In a previous version this was relative to input_dir, which is not right.
See [https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)](https://genome.sph.umich.edu/wiki/Regions_of_high_linkage_disequilibrium_(LD)) for a discussion.

* `output`: the base name of the output files. *This cannot be the same as the input!!!*

## 2 Overview of the workflow

The QC process consists of:

* removing duplicate markers;
* indentifying indviduals for whom there is discordant sex information;
* removing individuals with too high missingness or excessive heterozygosity;
* detecting whether there are any related individuals and removing enough to ensure that there are not related pairs;
* removing SNPs with too low MAF, or too high missingness, or anomalous HWE, or SNPs where there is a high differential missingness between cases and controls;
* a PCA of the resultant data is computed;
* a detailed report of the QC process is done.

## 3 Additional QC Parameters

The following parameters control QC

*  `sexinfo_available`: `true` or `false`. If we don't have sex information then we cannot do the check for discordant genotype. Note that it does not make sense (and is an error) to have sexinfo_available set to true when there is no X-chromosme data in the file;
*  `f_lo_male` and `f_hi_female`. Discordant sex genotype is done on the X-chromosome using the non-recombining parts. F, the in-breeding coefficient of the X-chromosome is computed. If F is above `f_lo_male`, the individual is predicted to be male; if F is below `f_hi_female`, the individual is predicted to be female. Anyone in between is flagged. These cut-off values are arbitrary and especially in large samples you are likely to find a range of F values. However, a large number of discrepant values probably indicates a sample mishandle error.  The PLINK default values (0.8 and 0.2) are the default parameters of the pipeline.
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


## 4 Performance parameters

There are three parameters that are important control performance. You probably won't need to change this, but feel free.

* `plink_mem_req` : specify in MB or GB how much memory your processes that use PLINK require;
*  `other_mem_req` : specify how much other processes need;
*  `max_plink_cores` : specify how many cores your PLINK processes can use. (This is only for those PLINK operations that are parallelisable. Some processes can't be parallelised and our workflow is designed so that for those processes only one core is used).


## 5 Output

A PDF report can be found in the output directory. This describes the process as well as what the inputs and outputs were.

Note that one issue that sometimes occurs in analysis is that there may over time be multple copies of the same file, perhaps with minor differences. To help version control, the PDF report captures the md5 checksums of inputs and outputs.
