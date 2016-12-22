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


