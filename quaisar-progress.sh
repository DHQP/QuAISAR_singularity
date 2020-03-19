#!/bin/sh -l

#$ -o quapro.out
#$ -e quapro.err
#$ -N quapro
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Wrapper for the main quaisar script to allow easy forwarding for log output and to attemp a progress bar
#
# Usage ./quaisar-pro.sh
#
# Output loction: default_config.sh_output_location/run_ID/
#
# Modules required: None
#
# v1.0 (03/17/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checking for proper number of arguments from command line
if [[ $# -lt 1  || $# -gt 9 ]]; then
	echo "If reads are in default location set in config file then"
  echo "Usage: ./quaisar_containerized.sh -c location_of_config_file -i location_of_reads 1|2|3|4 -o path_to_parent_output_folder_location name_of_output_folder [-a|r]"
	echo "filename postfix numbers are as follows 1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz"
  echo "You have used $# args"
  exit 3
fi
