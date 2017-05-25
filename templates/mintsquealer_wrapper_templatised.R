#!/usr/bin/env Rscript

# nextflow-template

# Wrapper script which calls mintsquealer function.
# mintsquealer function pre -processes and genotypes data from Illumina’s Infinium BeadChips (via crlmm) 

# Author: Lerato E. Magosi
# Contact: <magosil86 at gmail.com>
# License: (c) 12/12/16 Lerato E. Magosi, Released under GPL v.2

# Acknowledgements:
# crlmm developer: Matt Ritchie <mritchie at wehi.edu.au>

#Tip: make sure the samplesheet is in the same directory as the .idat files

# command line arguments
#args <- commandArgs(TRUE)

# mintsquealer run type (i.e. basic or advanced)
mintsquealer_run <- $mintsquealer_run_nf

# run on example dataset
run_demo_in <- $run_demo_in_nf


path2projdir_in <- $path2projdir_in_nf
log_filepath_in <- $log_filepath_in_nf
log_file_in <- $log_file_in_nf
path2idats_in <- $path2idats_in_nf
filename_samplesheet_in <- $filename_samplesheet_in_nf
cdf_name_in <- $cdf_name_in_nf 
run_quietly_in <- as.logical($run_quietly_in_nf)
show_advanced_params_in <- as.logical($show_advanced_params_in_nf)


if (mintsquealer_run == "advanced") {

    # set advanced params 
   num_markers_in <- as.numeric($num_markers_in_nf) 
   num_samples_in <- as.numeric($num_samples_in_nf)  
   set_snr_min_in <- as.numeric($set_snr_min_in_nf)
   set_call_method_in <- $set_call_method_in_nf 
   set_ids_in <- $set_ids_in_nf
   if (set_ids_in == "none") set_ids_in = NULL
   set_path_in <- $set_path_in_nf

   set_highDensity_in <- as.logical($set_highDensity_in_nf)
   set_sep_in <- $set_sep_in_nf 
   set_fileExt_in = list(green="Grn.idat", red="Red.idat")
   set_XY_in <- $set_XY_in_nf
   if (set_XY_in == "none") set_XY_in = NULL 
   set_trueCalls_in <- $set_trueCalls_in_nf
   if (set_trueCalls_in == "none") set_trueCalls_in = NULL
   set_copynumber_in <- as.logical($set_copynumber_in_nf)
   set_stripNorm_in <- as.logical($set_stripNorm_in_nf)
   set_useTarget_in <- as.logical($set_useTarget_in_nf)

 
   set_quantile_method_in <- $set_quantile_method_in_nf 
   set_mixtureSampleSize_in <- as.numeric($set_mixtureSampleSize_in_nf)
   set_fitMixture_in <- as.logical($set_fitMixture_in_nf)
   set_eps_in <- as.numeric($set_eps_in_nf)
   set_verbose_in <- as.logical($set_verbose_in_nf)
   set_sns_in <- $set_sns_in_nf
   if (set_sns_in == "none") set_sns_in = NULL
   set_probs_in = rep(1/3, 3)

   set_DF_in <- as.numeric($set_DF_in_nf)
   set_recallMin_in <- as.numeric($set_recallMin_in_nf) 
   set_recallRegMin_in <- as.numeric($set_recallRegMin_in_nf) 
   set_gender_in <- $set_gender_in_nf
   if (set_gender_in == "none") set_gender_in = NULL
   set_returnParams_in <- as.logical($set_returnParams_in_nf)
   set_badSNP_in <- as.numeric($set_badSNP_in_nf)

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
               num_markers = num_markers_in, num_samples = num_samples_in, set_snr_min = set_snr_min_in,
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
