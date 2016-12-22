#!/usr/bin/env Rscript


# Wrapper script which calls mintsquealer function.
# mintsquealer function pre -processes and genotypes data from Illumina’s Infinium BeadChips (via crlmm) 

# Author: Lerato E. Magosi
# Contact: <magosil86 at gmail.com>
# License: (c) 12/12/16 Lerato E. Magosi, Released under GPL v.2

# Acknowledgements:
# crlmm developer: Matt Ritchie <mritchie at wehi.edu.au>


# command line arguments
args <- commandArgs(TRUE)

# mintsquealer run type (i.e. basic or advanced)
mintsquealer_run <- args[1]

# run on example dataset
run_demo_in <- args[2]


path2projdir_in <- args[3]
log_filepath_in <- args[4]
log_file_in <- args[5]
path2idats_in <- args[6]
filename_samplesheet_in <- args[7]
cdf_name_in <- args[8]  
run_quietly_in <- args[9] 
show_advanced_params_in <- args[10]


if (mintsquealer_run == "advanced") {

    # set advanced params 
	num_markers_in <- args[11] 
	num_samples_in  <- args[12]  
	set_snr_min_in  <- args[13]
	set_call_method_in <- args[14] 
	set_ids_in  <- args[15] 
	set_path_in  <- args[16]

	set_highDensity_in <- args[17] 
	set_sep_in  <- args[18] 
	set_fileExt_in  = list(green = args[19], red = args[20])
	set_XY_in  <- args[21] 
	set_trueCalls_in  <- args[22] 
	set_copynumber_in  <- args[23]
	set_stripNorm_in  <- args[24] 
	set_useTarget_in  <- args[25]
 
	set_quantile_method_in  <- args[26] 
	set_mixtureSampleSize_in  <- args[27]
	set_fitMixture_in  <- args[28] 
	set_eps_in  <- args[29]
	set_verbose_in  <- args[30]
	set_sns_in  <- args[31] 
	set_probs_in = rep(1/3, 3)

	set_DF_in  <- args[32]
	set_recallMin_in  <- args[33] 
	set_recallRegMin_in  <- args[34] 
	set_gender_in  <- args[35]
	set_returnParams_in  <- args[36] 
	set_badSNP_in  <- args[37]

}


source("/Users/lmagosi/Downloads/crlmm_tutorials/mintsquealer/mintsquealer.R")



mintsquealer_wrapper <- function(mintsquealer_run) {

	if (mintsquealer_run == "basic") {

		mintsquealer(path2projdir = path2projdir_in, 
					 log_filepath = log_filepath_in, 
					 log_file = log_file_in,
					 path2idats = path2idats_in,
					 filename_samplesheet = filename_samplesheet_in,
					 cdf_name = cdf_name_in,  
					 run_quietly = run_quietly_in, 
					 show_advanced_params = show_advanced_params_in,
					 run_demo = run_demo_in)


	} else if (mintsquealer_run == "advanced") {

        mintsquealer(path2projdir = path2projdir_in, path2idats = path2idats_in, filename_samplesheet = filename_samplesheet_in, cdf_name = cdf_name_in, 
					num_markers = num_markers_in, num_samples = num_samples_in,  set_snr_min = set_snr_min_in,
					set_call_method = set_call_method_in, set_ids = set_ids_in, set_path = set_path_in,
					set_highDensity = set_highDensity_in, set_sep=set_sep_in, 
					set_fileExt = set_fileExt_in,
					set_XY = set_XY_in, set_trueCalls = set_trueCalls_in, set_copynumber = set_copynumber_in,
					set_stripNorm = set_stripNorm_in, set_useTarget = set_useTarget_in, 
					set_quantile.method = set_quantile_method_in, set_mixtureSampleSize = set_mixtureSampleSize_in,
					set_fitMixture = set_fitMixture_in, set_eps = set_eps_in, set_verbose = set_verbose_in, 
					set_sns = set_sns_in, set_probs = set_probs_in, set_DF = set_DF_in,
					set_recallMin = set_recallMin_in, set_recallRegMin = set_recallRegMin_in, set_gender = set_gender_in,
					set_returnParams = TRUE, set_badSNP = set_badSNP_in, show_advanced_params = show_advanced_params_in,
					run_demo = run_demo_in, run_quietly = run_quietly_in, log_filepath = log_filepath_in, log_file = log_file_in)

	}


}

# call the mintsquealer wrapper
mintsquealer_wrapper(mintsquealer_run)