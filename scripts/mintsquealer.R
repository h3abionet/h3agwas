# Goal: Pre -process and genotype data from Illumina’s Infinium BeadChips (via crlmm) 

# Author: Lerato E. Magosi
# Contact: <magosil86 at gmail.com>
# License: (c) 12/12/16 Lerato E. Magosi, Released under GPL v.2

# Acknowledgements:
# crlmm developer: Matt Ritchie <mritchie at wehi.edu.au>


mintsquealer_figlet <- function() {

	msg <- paste("~ ----------------------------------------------------------------------------------------------- ~",
				"|                     mintsquealer!, (C) 2016 Lerato E. Magosi, GPL v.2                            |",
				"|--------------------------------------------------------------------------------------------------|",
				"|                 ( Illumina{idat} --> genotype calls + confidence probs )                         |",
				"|--------------------------------------------------------------------------------------------------|",
				"|   A script that pre-processes and calls genotypes from Illumina's infinium bead chips via crlmm: |",
				"|                       https://github.com/magosil86/mintsquealer/                                 |",
				"~ ------------------------------------------------------------------------------------------------ ~",
	            sep = "\n")
	packageStartupMessage(msg)
}



# Description of mintsquealer params
# note: these are basically just params to crlmm's genotype.illumina function

# parameter                type               comment
# --------------------    -----------        ------------------------------

# path2projdir             character          set path to project directory
# path2idats               character          set path to .idat files
# filename_samplesheet     character          set path to Illumina sample sheet
# cdf_name                 character          set annotation package (crlmm has a valid list of annotation packages, see crlmm::validCdfNames())
# num_markers              numeric            set number of markers to process at a time
# num_samples              numeric            set number of samples to process at a time
# set_snr_min              numeric            set minimum signal-to-noise ratio used to filter out samples
# set_call_method          character          set calling algorithm e.g. crlmm or krlmm(only available for certain GWAS chips)
# set_ids                  character          set ids of probes to be read in (usually set to NULL)
# set_path                 character          set location of files to be read by genotype.illumina (usually just "." for current dir)
# set_highDensity          logical            used when sampleSheet is specified. If TRUE, array extensions ’\_A’, ’\_B’ in sampleSheet are replaced with ’R01C01’, ’R01C02 (usually set to FALSE)
# set_fileExt              list               set .idat file extensions for the Cy3 and Cy5 channels e.g. list(green="Grn.idat", red="Red.idat")
# set_XY                                      NChannelSet containing X and Y intensities (usually set to NULL)
# set_trueCalls            matrix             allows user to specify known Genotype calls(can contain some NAs) for a subset of samples and features (1 - AA, 2 - AB, 3 - BB)
# set_copynumber           logical            Whether to store copy number intensities with SNP output.
# set_stripNorm            logical            Whether the data should be strip-level normalized?
# set_useTarget            logical            (only used when stripNorm=TRUE). Should the reference HapMap in- tensities be used in strip-level normalization?
# set_quantile.method      character          set quantile normalization method to use (’within’ or ’between’ channels)
# set_mixtureSampleSize                       Sample size to be use when fitting the mixture model.
# set_fitMixture           logical            Whether to fit per-array mixture model
# set_eps                                     Stop criteria
# set_verbose              logical            Whether to print descriptive messages during processing.
# sns                                         set sample identifiers. If missing, the default sample names are basename(filenames)
# set_probs                numeric            vector with priors for AA, AB and BB
# set_DF                                      number of degrees of freedom to use with t-distribution
# set_recallMin            integer            Minimum number of samples for recalibration.
# set_recallRegMin         integer            Minimum number of SNP’s for regression          
# set_gender               integer vector     ( male = 1, female = 2 ) or missing, with same length as filenames. If missing, the gender is imputed.
# set_returnParams         logical            Return recalibrated parameters from crlmm
# set_badSNP               numeric            Threshold to flag bad SNP (affects batchQC)
# show_advanced_params     logical            display values assigned to params
# run_demo                 logical            run on example dataset
# run_quietly              logical            Whether to run mintsquealer verbosely or quietly
# log_filepath             character          path2logfile
# log_file                 character          log filename

# ------------------------------------------------------------------------------------------



# basic usage:

# mintsquealer(path2projdir = "/Users/lmagosi/Downloads/crlmm_tutorials/", 
#              log_filepath = "/Users/lmagosi/Downloads/crlmm_tutorials/", 
#              log_file = "mintsquealer_log",
#              path2idats = base::system.file("idatFiles", package = "hapmap370k"),
#              filename_samplesheet = "samples370k.csv",
#              cdf_name = "human370v1c",  
#              run_quietly = F, 
#              show_advanced_params = TRUE)


# Example analysis:
# mintsquealer(path2projdir = "/Users/lmagosi/Downloads/crlmm_tutorials/", 
#              log_filepath = "/Users/lmagosi/Downloads/crlmm_tutorials/", 
#              log_file = "mintsquealer_log",
#              show_advanced_params = TRUE,
#              run_demo = T) 



mintsquealer <- function(path2projdir, path2idats, filename_samplesheet, cdf_name, 
                            num_markers = 150e3, num_samples = 500,  set_snr_min = 5,
                            set_call_method = "crlmm", set_ids = NULL, set_path = ".",
                            set_highDensity = FALSE, set_sep="_", 
                            set_fileExt = list(green="Grn.idat", red="Red.idat"),
                            set_XY = NULL, set_trueCalls = NULL, set_copynumber = TRUE,
                            set_stripNorm = TRUE, set_useTarget = TRUE, 
                            set_quantile.method = "between", set_mixtureSampleSize = 10^5,
                            set_fitMixture = TRUE, set_eps = 0.1, set_verbose = TRUE, 
                            set_sns = NULL, set_probs = rep(1/3, 3), set_DF = 6,
                            set_recallMin = 10, set_recallRegMin = 1000, set_gender = NULL,
                            set_returnParams = TRUE, set_badSNP = 0.7, show_advanced_params = FALSE,
                            run_demo = FALSE, run_quietly = TRUE, log_filepath, log_file) {


	# generate timestamp for logfile
	t_logfile <- base::Sys.time()
	timeStamp_logfile <- base::strftime(t_logfile,"%Y_%m_%d__%H_%M_%S")

	# send console output to logfile
	connection <- base::file(base::paste0(log_filepath, "/", log_file, "_", timeStamp_logfile, ".Rlog"))
	base::sink(connection, append=TRUE)
	base::sink(connection, append=TRUE, type="message")

   # ----------------------------------------------------------------------------

    # call figlet function
    mintsquealer_figlet()
    base::writeLines("")

   # ----------------------------------------------------------------------------

	if (run_demo) {

		# set path to illumina raw data (.idat) files
		path2idats <- base::system.file("idatFiles", package = "hapmap370k")

		# Illumina sample sheets contain information for reading raw IDAT files
		# See BeadStudio Genotyping guide, Appendix A for mor info
		filename_samplesheet <- "samples370k.csv"

		# set annotation package
		cdf_name <- "human370v1c"

	}    

   # ----------------------------------------------------------------------------
	

	
	# Parameter_values
	base::writeLines("Basic parameter values for current mintsquealer run: ")
	base::writeLines("")
	
	base::print(base::paste0("path2projdir (path to project directory): ", path2projdir))
    base::writeLines("")

	base::print(base::paste0("path2idats (path to illumina raw data {.idat} files): ", path2idats))
    base::writeLines("")
	
	base::print(base::paste0("filename_samplesheet (Illumina sample sheet): ", filename_samplesheet))
    base::writeLines("")

	base::print(base::paste0("cdf_name (annotation package) : ", cdf_name))
    base::writeLines("")
    
    base::print(base::paste0("num_markers (number of markers to process at a time): ", num_markers))
    base::writeLines("")
    
    base::print(base::paste0("num_samples (number of samples to process at a time): ", num_samples))
    base::writeLines("")

	base::print(base::paste0("set_snr_min (minimum signal-to-noise ratio threshold to filter out samples): ", set_snr_min))
    base::writeLines("")


   # ----------------------------------------------------------------------------
    
    if (show_advanced_params) {

	base::writeLines("Advanced parameter values for current mintsquealer run: ")
	base::writeLines("")

	base::writeLines("For a description of the advanced parameters: ?genotype.Illumina ")
	base::writeLines("")

	base::print(base::paste0("ids: ", set_ids))
    base::writeLines("")

	base::print(base::paste0("path: ", set_path))
    base::writeLines("")

	base::print(base::paste0("highDensity: ", set_highDensity))
    base::writeLines("")

	base::print(base::paste0("sep: ", set_sep))
    base::writeLines("")

	base::print(base::paste0("fileExt: ", set_fileExt))
    base::writeLines("")

	base::print(base::paste0("XY: ", set_XY))
    base::writeLines("")

	base::print(base::paste0("trueCalls: ", set_trueCalls))
    base::writeLines("")

	base::print(base::paste0("copynumber: ", set_copynumber))
    base::writeLines("")

	base::print(base::paste0("stripNorm: ", set_stripNorm))
    base::writeLines("")

	base::print(base::paste0("useTarget: ", set_useTarget))
    base::writeLines("")

	base::print(base::paste0("quantile.method: ", set_quantile.method))
    base::writeLines("")
 
 	base::print(base::paste0("mixtureSampleSize: ", set_mixtureSampleSize))
    base::writeLines("")

	base::print(base::paste0("fitMixture: ", set_fitMixture))
    base::writeLines("")

	base::print(base::paste0("eps: ", set_eps))
    base::writeLines("")

	base::print(base::paste0("verbose: ", set_verbose))
    base::writeLines("")

	base::print("seed: is randomly generated and reported at the end of the run")
    base::writeLines("")

	base::print(base::paste0("set_sns: ", set_sns))
    base::writeLines("")

	base::print(base::paste0("probs: ", set_probs))
    base::writeLines("")

	base::print(base::paste0("DF: ", set_DF))
    base::writeLines("")

	base::print(base::paste0("recallMin: ", set_recallMin))
    base::writeLines("")

	base::print(base::paste0("recallRegMin: ", set_recallRegMin))
    base::writeLines("")

	base::print(base::paste0("gender: ", set_gender))
    base::writeLines("")

	base::print(base::paste0("returnParams: ", set_returnParams))
    base::writeLines("")

	base::print(base::paste0("badSNP: ", set_badSNP))
    base::writeLines("")
        
   
    }


   # ----------------------------------------------------------------------------


	# Load libraries ----------------

	suppressPackageStartupMessages(library(ff))
	options(ffcaching="ffeachflush") # helps to prevent system stalls on large writes
	suppressPackageStartupMessages(library(crlmm))
	suppressPackageStartupMessages(library(Biobase))

	# example dataset
	if (run_demo) library(hapmap370k) # example dataset: 40 HapMap samples analyzed using Illumina’s 370k Duo BeadChips
	if (run_demo) library(human370v1cCrlmm) # chip-specific model parameters and basic SNP annotation information for hapmap370k
    
	# Set project dir ----------------

	# create project dir. to store results

	t <- base::Sys.time()
	timeStamp <-  base::strftime(t,"%Y_%m_%d__%H_%M_%S")
	outdir <- base::paste0(path2projdir, "crlmm_project_", timeStamp)
	ff_files <- base::paste0(base::file.path(outdir, "ff_files"))
	base::dir.create(ff_files, recursive=TRUE)

	# Supplying the outdir to ldpath ensures that all ff files generated during genotyping
	# are stored in the outdir
	oligoClasses::ldPath(ff_files)


	# Set genotyping environment ----------------

	# set number of markers to process at a time
	oligoClasses::ocProbesets(num_markers)

	# set number of samples to process at a time
	oligoClasses::ocSamples(num_samples)


	# Load data ----------------


	# load sample sheet
	# Illumina sample sheets contain information for reading raw IDAT files
	# See BeadStudio Genotyping guide, Appendix A for mor info

	samplesheet <- utils::read.csv(base::file.path(path2idats, filename_samplesheet), header = TRUE, as.is = TRUE)

	# inspect sample sheet	 
	if (run_quietly == FALSE) { base::writeLines("head(samplesheet): "); base::print(utils::head(samplesheet)); base::writeLines("")}

	# set array names i.e. absolute path to the idats filenames without the file extension
	array_names <- base::file.path(path2idats, base::unique(samplesheet[, "SentrixPosition"]))

	# check if all of the green and red IDAT files exist
	# check green intensity files
	if(!(base::all(base::file.exists(base::paste0(array_names, "_Grn.idat"))))) base::stop("Check whether all the green IDAT files exist: '*_Grn.idat'") 

	# check red intensity files
	if(!(base::all(base::file.exists(base::paste0(array_names, "_Red.idat"))))) base::stop("Check whether all the red IDAT files exist: '*_Red.idat'")

	# set array info i.e. barcode and position
	array_barcode <- NULL
	array_info <- base::list(barcode = array_barcode, position = "SentrixPosition")
	if (run_quietly == FALSE) { base::writeLines("array_info: "); base::print(array_info)}


	# Extract Red and Green intensities from Illumina binary (.idat) files ----------------


	# create an instance/object of the NChannelSet class to store Red and Green intensities
	base::system.time(rg <- crlmm::readIdatFiles(samplesheet, path = path2idats, arrayInfoColNames = array_info, saveDate = TRUE))

	# check class and slot names of the rg object
	if (run_quietly == FALSE) { 
	    base::writeLines("class and slot names of the rg object: ")
	    base::print(base::class(rg))
	    base::writeLines("")
	    base::print(methods::slotNames(rg))
	    base::writeLines("")
	    }

	# check dimensions and structure of the rg object
	if (run_quietly == FALSE) {base::writeLines("dim(rg): "); base::print(base::dim(rg))}
	base::writeLines("")
	base::writeLines("Finished extracting Red and Green intensities from .idat files.")
	base::writeLines("The Structure of the NChannelSet class object where the Red and Green intensities are stored is shown below: ")
	base::writeLines("")
	utils::str(rg)

	# inspect red intensities for the first five samples
	# note: exprs function accesses assay data for objects derived from the 'eSet-class' (found in Biobase)
	if (run_quietly == FALSE) {
	    base::writeLines("")
	    base::writeLines("channelNames(rg): ")
	    base::print(Biobase::channelNames(rg))
	    base::writeLines("")
	    base::writeLines("inspect red intensities for the first five samples: ")
	    base::print(Biobase::exprs(Biobase::channel(rg, "R"))[1:5, 1:5])

	# inspect green intensities for the first five samples
	    base::writeLines("")
	    base::writeLines("inspect green intensities for the first five samples: ")
	    base::print(Biobase::exprs(Biobase::channel(rg, "G"))[1:5, 1:5])
	    base::writeLines("")
	    
	    }

	# Retrieve information on experimental phenotypes
	# note: 'pData' returns a data.frame with samples as rows, variables as columns
	pd <- Biobase::pData(rg)
	if (run_quietly == FALSE) {base::writeLines("Retrieve information on experimental phenotypes: "); base::print(utils::head(pd)); base::writeLines("")}



	# Determine number of array batches ----------------


	# get scanning dates for each array
	rg@protocolData
	scandatetime <- base::strptime(Biobase::protocolData(rg)[["ScanDate"]], "%m/%d/%Y %H:%M:%S %p")
	datescanned <- base::substr(scandatetime, 1, 10)
	if (run_quietly == FALSE) {base::writeLines("get scanning dates for each array: "); base::print(utils::head(datescanned))}

	# convert datescaneed to a factor
	scanbatch <- base::factor(datescanned)

	# replace levels of factor variable datescanned with an integer sequence
	base::levels(scanbatch)
	base::writeLines("")
	base::levels(scanbatch) <- base::seq_along(base::levels(scanbatch))
	base::levels(scanbatch)

	# convert scanbatch to a vector
	scanbatch <- base::as.numeric(scanbatch)
	if (run_quietly == FALSE) {base::writeLines("convert scanbatch to a vector: "); base::print(scanbatch)}


	# Check for arrays with poor signal ----------------


	# create a plots dir
	plots_dir <- base::paste0(base::file.path(outdir, "plots"))
	base::dir.create(plots_dir, recursive=TRUE)

	# generate plot showing array signal intensities for each sample

	# get chip name
	name_chip <- base::substr(filename_samplesheet, 1, base::nchar(filename_samplesheet) - 4)

	# set filename for plot
	filename_signal_qlty_arrays <- base::file.path(plots_dir, base::paste0("graph_array_signal_qlty_", name_chip, "_", timeStamp, ".tif"))

	# get number of samples
	names_samples <- base::seq(base::dim(rg)[["Samples"]])

	# plot
	grDevices::tiff(filename_signal_qlty_arrays , width = 17.35, height = 23.35, units = "cm", res = 300, compression = "lzw", pointsize = 14, colortype = "true", family = "Helvetica")
	graphics::par(mfrow = base::c(2, 1), mai = base::c(0.4, 0.4, 0.4, 0.1), oma = base::c(1, 1, 0, 0))
	graphics::boxplot(base::log2(Biobase::exprs(Biobase::channel(rg, "R"))), xlab = "Array", ylab = "", names = names_samples, main = "Red channel", outline = FALSE, las = 2)
	graphics::boxplot(base::log2(Biobase::exprs(Biobase::channel(rg, "G"))), xlab = "Array", ylab = "", names = names_samples, main = "Green channel", outline = FALSE, las = 2)
	graphics::mtext(base::expression(log[2](intensity)), side = 2, outer = TRUE)
	grDevices::dev.off()


	# Pre -processing and Genotyping ----------------
	base::writeLines("")
    base::writeLines("Pre -processing and Genotyping via CRLMM: ")
    base::writeLines("")

	# set annotation package

	# check validity of the submitted annotation package
	cdf_name %in% crlmm::validCdfNames()
	if (!(cdf_name %in% crlmm::validCdfNames())) { base::print(base::paste0(cdf_name, " is an invalid annotation package for crlmm")); base::writeLines(""); base::writeLines("Here is a list of valid annotation packages: "); base::print(crlmm::validCdfNames()); base::stop(" consider one of the valid annotation packages listed above. ")}

	# create a character vector that specifies the batch for each array. 
	char_scanbatch <- base::as.character(scanbatch)

	# use random number generator to select a seed value
	vec_seed_numbers <- 1:4e06
	pseudo_random_seed <- base::sample(vec_seed_numbers, 1)
	if (run_quietly == FALSE) {writeLines("pseudo_random_seed: "); base::print(pseudo_random_seed); base::writeLines("")}

	# genotype.Illumina is a wrapper that calls the functions: preprocessInf and genotypeInf
	# preprocessInf, reads and normalizes intensities from .idat files and returns ff objects (contain params for the mixture model)
	# genotypeInf, does the genotyping via the crlmm or krlmm algorithms


	# call genotypes
	cnSet <- crlmm::genotype.Illumina(sampleSheet = samplesheet,
								arrayNames = array_names,
								arrayInfoColNames = array_info,
								cdfName = cdf_name,
								batch = char_scanbatch,
								seed = pseudo_random_seed,
								call.method = set_call_method,
								SNRMin = set_snr_min,
								saveDate = TRUE,
								
								# advanced params
								ids = set_ids, path = set_path, highDensity = set_highDensity,
								sep = set_sep, fileExt = set_fileExt, XY = set_XY, trueCalls = set_trueCalls,
								copynumber = set_copynumber, stripNorm = set_stripNorm, 
								useTarget = set_useTarget, quantile.method = set_quantile.method,
								mixtureSampleSize = set_mixtureSampleSize, fitMixture = set_fitMixture,
								eps = set_eps, sns = set_sns, verbose = set_verbose, probs = set_probs,
								DF = set_DF, recallMin = set_recallMin, recallRegMin = set_recallRegMin,
								gender = set_gender, returnParams = set_returnParams, badSNP = set_badSNP)



	# check if the cnSet object is an instance of the CNSet class
	if (run_quietly == FALSE) {base::writeLines("class(cnSet): "); base::print(base::class(cnSet))}

	# inspect the structure of the cnSet object
	base::writeLines("Finished pre -processing and Genotyping via CRLMM.")
	base::writeLines("The Structure of the CNSet class object where the genotypes are stored is shown below: ")
	base::writeLines("")
	utils::str(cnSet)

	# check dimensions and slot names of the cnSet object
	if (run_quietly == FALSE) {
	    base::writeLines("")
	    base::writeLines("check dimensions and slot names of the cnSet object: ")
	    base::print(base::dim(cnSet))
	    base::writeLines("")
	    base::writeLines("slotNames(cnSet): ")
	    base::print(methods::slotNames(cnSet))
	    

	# inspect genotype calls in first five samples
	    base::writeLines("")
	    base::writeLines("inspect genotype calls in first five samples: ")
	    base::print(oligoClasses::calls(cnSet)[1:10, 1:5])

	# inspect confidence probabilities in first five samples
	    base::writeLines("")
	    base::writeLines("confidence probabilities in first five samples: ")
	    base::print(oligoClasses::i2p(oligoClasses::confs(cnSet)[1:10, 1:5]))
	    base::writeLines("")
	}


	# generate plot showing signal-to-noise ratio per array to check for batch effects
	# note: different symbols are used for arrays scanned on different days (i.e. pch is set to scanbatch)

	# inspect the structure of the SNR in the cnSet object
	base::writeLines("inspect the structure of the SNR in the cnSet object: ")
	utils::str(cnSet[["SNR"]])
	base::writeLines("")

	# set filename for plot
	filename_snr_per_array <- base::file.path(plots_dir, base::paste0("graph_snr_per_array", name_chip, "_", timeStamp, ".tif"))


	invisible(oligoClasses::open(cnSet$SNR))
	snr <- cnSet$SNR[]
	oligoClasses::close(cnSet$SNR)

	# plot
	grDevices::tiff(filename_snr_per_array, width = 17.35, height = 23.35, units = "cm", res = 300, compression = "lzw", pointsize = 14, colortype = "true", family = "Helvetica")                
	graphics::plot(snr, pch = scanbatch, xlab = "Array", ylab = "SNR", main = "Signal-to-noise ratio per array", las = 2)
	grDevices::dev.off()

	# generate histogram showing signal-to-noise ratio
	# note: SNR measures array quality (i.e. separation of diallelic genotypic clusters at polymorphic loci)
	#       Small SNR values can indicate possible problems with the DNA.
	#       SNR below 5 may indicate poor sample quality for Affymetrix data.
	#       SNR below 25 may indicate poor sample quality for Illumina data.


	# set filename for plot
	filename_snr <- base::file.path(plots_dir, base::paste0("graph_snr_", name_chip, "_", timeStamp, ".tif"))

	# plot
	grDevices::tiff(filename_snr, width = 17.35, height = 23.35, units = "cm", res = 300, compression = "lzw", pointsize = 14, colortype = "true", family = "Helvetica")                
	graphics::hist(snr, breaks = 25, xlim = range(snr), xlab = "SNR", main = "Signal-to-noise ratio", density = TRUE)
	grDevices::dev.off()


	# Write out calls and confidence scores  ----------------


	results <- base::paste0(base::file.path(outdir, "results"))
	base::dir.create(results, recursive=TRUE)

    calls_file <- base::file.path(results, "crlmm_calls.txt")
	calls_ffdf <- ff::as.ffdf(oligoClasses::calls(cnSet))
	write.table.ffdf(calls_ffdf, file = calls_file)
	oligoClasses::close(oligoClasses::calls(cnSet))


	#library(MASS)
	confs_file <- base::file.path(results, "crlmm_confs.txt")
	confs_matrix <- base::t(base::as.matrix(oligoClasses::i2p(oligoClasses::confs(cnSet)[,]), rownames.force = TRUE))
	utils::write.table(confs_matrix, file = confs_file)
	oligoClasses::close(oligoClasses::confs(cnSet))

	

	# Record session info  ----------------


	# get session info for record keeping
	base::writeLines("")
	base::writeLines("session info: ")
	base::writeLines("")
	base::print(utils::sessionInfo())


	# And that's all folks!

	# restore output to console
	base::sink() 
	base::sink(type="message")


    }

