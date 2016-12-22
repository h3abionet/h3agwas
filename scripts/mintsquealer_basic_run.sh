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
MINTSQUEALER_RUN="basic"

# run demo analysis
RUN_DEMO_IN=TRUE


# set basic params
PATH2PROJDIR_IN="/Users/lmagosi/Downloads/crlmm_tutorials/"
LOG_FILEPATH_IN="/Users/lmagosi/Downloads/crlmm_tutorials/"
LOG_FILE_IN="mintsquealer_log"
PATH2IDATS_IN="base::system.file('idatFiles', package = 'hapmap370k')"
FILENAME_SAMPLESHEET_IN="samples370k.csv"
CDF_NAME_IN="human370v1c"  
RUN_QUIETLY_IN=FALSE 
SHOW_ADVANCED_PARAMS_IN=TRUE


# set advanced params 
NUM_MARKERS_IN=150e3 
NUM_SAMPLES_IN=500  
SET_SNR_MIN_IN=5
SET_CALL_METHOD_IN="crlmm" 
SET_IDS_IN=NULL 
SET_PATH_IN="."

# SET_HIGHDENSITY_IN=FALSE 
# SET_SEP_IN="_"
# GREEN="Grn.idat" 
# RED="Red.idat"
# SET_XY_IN=NULL 
# SET_TRUECALLS_IN=NULL 
# SET_COPYNUMBER_IN=TRUE
# SET_STRIPNORM_IN=TRUE 
# SET_USETARGET_IN=TRUE

# SET_QUANTILE_METHOD_IN="between" 
# SET_MIXTURESAMPLESIZE_IN=10^5
# SET_FITMIXTURE_IN=TRUE 
# SET_EPS_IN=0.1
# SET_VERBOSE_IN=TRUE
# SET_SNS_IN=NULL 


# SET_DF_IN=6
# SET_RECALLMIN_IN=10 
# SET_RECALLREGMIN_IN=1000 
# SET_GENDER_IN=NULL
# SET_RETURNPARAMS_IN=TRUE
# SET_BADSNP_IN=0.7




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
					$SHOW_ADVANCED_PARAMS_IN << EOD


EOD

echo *******************************************
echo "Finished at: "`date`
echo *******************************************
exit 0


