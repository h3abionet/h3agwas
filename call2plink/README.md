
# 7. Converting from Illumina genotyping reports in TOP/BOTTOM or Forward foramtformat

This workflow is run by 

```nextflow run call2plink```

and converts from an Illumina TOP/BOTTOM or FORWARD call file. Together with auxiliary input data, this file is first converted into a raw PLINK file and then the PLINK file is aligned to a strand, and then convered into binary PLINK format. This process can take a very long time.

If you don't understand these formats, the bad news is that you really should. See  S. C. Nelson, K. F. Doheny, C. C. Laurie, and D. B. Mirel, "Is 'forward' the same as 'plus'?...and other adventures in SNP allele nomenclature." Trends in Genetics, vol. 28, no. 8, pp. 361-363, 2012. http://dx.doi.org/10.1016/j.tig.2012.05.002  The good news is that you are no talone

Essentially the problem since the reference genome changes over time, what is on the forward strand of one reference could become the reverse in the next. Not likely but could and does happen. Thus what someone sees as a SNP with A/C alleles could become a SNP with T/G alleles etc. For SNPs with A/C alleles we can easily see when something's flipped but if the allele is an A/T or C/G allele we can't differentiate between alternate alleles and reverse complements.

Two common methods that are used to disambiguuate this are to call
* with respect to the entry in dbSNP -- usually the way in which the discoverer reported it. This submission to dbSNP will contain the flanking regions of the SNP. Usually this is will be in the smae orientation as the reference genome, but often is not 
* Illumina's TOP/BOTTOM format which uses the SNP and/or flanking region (see the reference given)

 
This process is expensive because:

* the top/bottom or forward file is a very bulky and inefficient format
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

* `manifest`: The chip manifest file. This is crucial. You can find examples here: https://emea.illumina.com/products/by-type/microarray-kits.html, but you may have to ask Illumina for the correct vrsion.

* `chipdescription`: this is a standard file produced by Illumina for your chip which contains not only the chromosome/coordinate of each SNP but also the genomic position (measured in centimorgans). If you don't have this -- give the manifest file. All will work except your bim files will not contain genonomic positoin

* `samplesheet`: This is Excel spreadsheet or CSV (comma-separated only) that Illumina or a genotyping centre provides which details each perfson in the study for whom you have genotyping results. If you don't have it, you can set this variable to 0 or the empty string, in which case the output PLINK fam file will have unknown values for sex and phenotype.  Alternatively, if you don't have it, you can make your own. Note that if the suffix is ".xls" or ".xlsx" the code assumes it's an Excel file otherwise a comma-separated set of values.

There are three columns that are important: "Institute Sample Label", "Manifest Sex" and "Batch Comment". These must be there. The _label_ is the ID of the person. In the current workflow this ID is used for both the FID and IID. If you have a family study you may need to manually change the fam file.


Please note that *we expect all entries in the sample IDs etc to be alphanumeric 0-9, Latin letters (NO accents!), underscore, space, hyphen*. The code may break otherwise.

* `idpat`. Default is "0" or "" (ignore). By default, we use the sample ID as found in the genotype report and sample sheet. PLINK fam files require a double barrelled name (FID IID) -- we just double the ID as found. However, this may not be ideal since
the Illumina IDs in the sample ID are typically a long string some of  the components of which will not be useful when you are analysing the result. You can change the sample ID by providing a Python-style regular expression which decribes the components. The regex groups describe the components. If there is one group, it is doubled. If there are two groups, then those become the FID and IID. Only one or two groups are permissible. 

For example, suppose the ID as found in the Illumina input data is `WG0680781-DNA_A02_ABCDE`, if you use ".*_(.+)" as the idpat, then the FID IID used would be ABCDE ABCDE. If you used "(\\w+)_DNA_(\\w+)_" then the FID IIS used would be "WG0680781 A02". Note how we need to escape the backslash twice.


Unfortunately we experience that genotyping centres have different formats and that you can even get the same centre changing the labels of columns of the report. Using the `sheet_columns` parameter you can make adjustmens.

* `params.sheet_columns`: this should be a file name. The file should explain what the column labels in your sample sheet are. The format is shown in the example below, where the default values are given (if you are happy with all of them you don't need the `sheet_columns` parameter -- if you are happy with some of them only put the ones you want to change). Here we are saying that the _sex_ as provided by the manifest is found in a column called "Manifest Sex", the sample is found in a column "Institute Sample Label" and so on. The first four are required by the workflow. If you don't have batch information, you can define `batch` as 0

````
sex=Manifest Sex
sample_label=Institute Sample Label
plate=Sample Plate
well=Well
batch=Batch Comment
````

* `output_align`. This can be one of three values: _dbsnp_, _ref_, and _db2ref_. dnsnp and ref assume that the input is in TOP/BOT format. If dbsnp, the output will be aligned to the dbSNP report, if "ref", the output will be aligned to a given reference strand. Many of the SNPs will be flipped (e.g. an A/C SNP will become G/T; and A/T SNP will become T/A).   _db2ref_ assumes the input is in FORWARD format and aligns to to the given reference genome.

* `strandreport`: This is an Illumina-style strand report. It is not needed if you choose "ref" above, but it is needed for the others.

* `refererence`: This is the name of a file that gives the reference allele for each SNP on the chip.  This is only useful if the "ref" option is used for `output_align`, and is optional in this case. Note that the difference between aligning and the use of the reference. Aligning will decide which strand of the reference genome as found in the Illumina genotyping teh alleles can be found on. For example, if the genotyping report gives the two options as A and C, aligning checks whether this is A and C on the + strand of the reference genome or (and so will be A and C in the output bim file) or whther this is A and C on the $-$ strand of the reference genome and so should be reported as T and G. This is done using information in the chip manifest file. The second step is to know which allele is the reference allele and which is the alternate allele.

A reference file suitable for the H3A chip can be found here http://www.bioinf.wits.ac.za/data/h3agwas/. Two formats for the reference file are supported: (1) simple -- two columns, no header row, the first column ins the SNP ID, the second the reference allele; and (2) complex -- >= two columns, one header row, the header row must contain a column label SNP and a column label Base for the SNP ID and reference allele respectively, all other columns are ignored.

* `batch_col`: For this workflow, the `batch_col` parameter is a column in the `samplesheet` that should be used to extract out out the 6-th column of the `fam` file, or the phenotype. This allows you do do batch analysis. Of course, you can choose anything you like to be the "batch". The default value is 0, which means just set the 6-th column of the fam file to -9.  One special case: If the contents of the column is of the form "Batch n", then only the _n_ is returned.

* `samplesize`: This was  included mainly for development purposes but _perhaps_ might be helpful to some users. This allows you sample only the first _n_ people in each genotype report. This allows you to extract out a small subset of the data for testing purposes. The default is 0, which means that *all* individuals will be generated.

* `indiv_memory_req`: default is 2GB -- good enough for 100 samples, 2.3 million SNPs per files (you may need to change this if you have one massive file rather than many big files)
* `combined_memory_req` : default is 8GB -- probably good enough for 15000 samples, 2.3millions SNPs in total

* `time_req`: default is 12h -- as above

### Advanced features for sample handling

* `mask`: This is a file of sample IDs that you want excluded from your data. These should be IDs given in the _Institute Sample Label_ of the sample sheet. The file should contain at least one column, possibly with other columns white-space delimited. Only the first column is used the other columns are ignored.

* `replicates`: This is a file of sample IDs that are biological replicates. You will often include biological replicates for genotyping -- the label as given in the _Institute Sample Label_ column will of course be different, but once you have extracted out the sample ID using the _idpat_ field above, all the replicates for the same individual will then have the same sample id. For samples that have replicates you should choose one of the samples to be the canonical one and then identify the others as being the replicates with the labels 

* `newpat`: This is experimental and only should be used with care. Suppose, completely hypothetically, there's a sample mix-up. The person you called "X3RTY" is actually "UYT0AV" who is actually "R2D2" and so on. You can fix your sample-sheet but the genotyping calls still have the same (wrong) values. If you have chosen your _Institute Sample_Label_ so that it contains both the ID and the plate and well then our scripts can help you. If not, good luck to you.
    * set _idpat_ to a regular expression that gives the plate ID, well as the FID and IID. This will give you the id uniquely determined by the plate and the well.
    * fix your sample sheet -- just fix the _Institute Sample Label_ field. 
    * Choose your _newpat_ as a regular expression that extracts out the correct ID from this. Our workflow will use the plate and well in your sample sheet to produce to match the plate/well from the genotype calling phase to the correct ID. (If you look in the working directory of the fixFam process the .command.out file will show you all the matches)


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
