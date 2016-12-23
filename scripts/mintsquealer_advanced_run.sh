#!/bin/bash

# Goal: Pre -processes and genotypes data from Illumina’s Infinium BeadChips (via crlmm)


# qsub options

# set environment

# set job name

# set project and queue


# user definition on SGE

echo *******************************************
echo "Goal: Pre -processes and genotypes data from Illumina’s Infinium BeadChips (via crlmm)"
echo -------------------------------------------
echo "Run on host: "`hostname`
echo "Operating system: "`uname -s`
echo "Username: "`whoami`
echo "Started at: "`date`
echo *******************************************


# set run type
MINTSQUEALER_RUN="advanced"

# run demo analysis
RUN_DEMO_IN=TRUE


# set basic params
PATH2PROJDIR_IN="/Users/lmagosi/Downloads/crlmm_tutorials/"
LOG_FILEPATH_IN="/Users/lmagosi/Downloads/crlmm_tutorials/"
LOG_FILE_IN="mintsquealer_log"
PATH2IDATS_IN="/Library/Frameworks/R.framework/Versions/3.1/Resources/library/hapmap370k/idatFiles"
FILENAME_SAMPLESHEET_IN="samples370k.csv"
CDF_NAME_IN="human370v1c"  
RUN_QUIETLY_IN=FALSE 
SHOW_ADVANCED_PARAMS_IN=TRUE


# set advanced params 
NUM_MARKERS_IN=150000
NUM_SAMPLES_IN=500  
SET_SNR_MIN_IN=5
SET_CALL_METHOD_IN="crlmm" 
SET_IDS_IN="none" 
SET_PATH_IN="."

SET_HIGHDENSITY_IN=FALSE 
SET_SEP_IN="_"
SET_XY_IN="none" 
SET_TRUECALLS_IN="none" 
SET_COPYNUMBER_IN=TRUE
SET_STRIPNORM_IN=TRUE 
SET_USETARGET_IN=TRUE

SET_QUANTILE_METHOD_IN="between" 
SET_MIXTURESAMPLESIZE_IN=100000
SET_FITMIXTURE_IN=TRUE 
SET_EPS_IN=0.1
SET_VERBOSE_IN=TRUE
SET_SNS_IN="none"


SET_DF_IN=6
SET_RECALLMIN_IN=10 
SET_RECALLREGMIN_IN=1000 
SET_GENDER_IN="none"
SET_RETURNPARAMS_IN=TRUE
SET_BADSNP_IN=0.7




# set paths
SCRIPTS_DIR="/Users/lmagosi/Downloads/crlmm_tutorials/mintsquealer"

cd $SCRIPTS_DIR

Rscript --vanilla $SCRIPTS_DIR/mintsquealer_wrapper.R $MINTSQUEALER_RUN \
                    $RUN_DEMO_IN \
					$PATH2PROJDIR_IN \
					$LOG_FILEPATH_IN \
					$LOG_FILE_IN \
					$PATH2IDATS_IN \
					$FILENAME_SAMPLESHEET_IN \
					$CDF_NAME_IN \
					$RUN_QUIETLY_IN \
					$SHOW_ADVANCED_PARAMS_IN \
					$NUM_MARKERS_IN \
					$NUM_SAMPLES_IN \
					$SET_SNR_MIN_IN \
					$SET_CALL_METHOD_IN \
					$SET_IDS_IN \
					$SET_PATH_IN \
					$SET_HIGHDENSITY_IN \
					$SET_SEP_IN \
					$SET_XY_IN \
					$SET_TRUECALLS_IN \
					$SET_COPYNUMBER_IN \
					$SET_STRIPNORM_IN \
					$SET_USETARGET_IN \
					$SET_QUANTILE_METHOD_IN \
					$SET_MIXTURESAMPLESIZE_IN \
					$SET_FITMIXTURE_IN \
					$SET_EPS_IN \
					$SET_VERBOSE_IN \
					$SET_SNS_IN \
					$SET_DF_IN \
					$SET_RECALLMIN_IN \
					$SET_RECALLREGMIN_IN \
					$SET_GENDER_IN \
					$SET_RETURNPARAMS_IN \
					$SET_BADSNP_IN << EOD


EOD

echo *******************************************
echo "Finished at: "`date`
echo *******************************************
exit 0


