#!/bin/sh -l

#$ -o run_SPAdes.out
#$ -e run_SPAdes.err
#$ -N run_SPAdes
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Runs SCCmesFinder on isolate to find mec types in staph samples
#
# Usage: ./run_SCCmecFinder.sh sample_name	run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/SCCmec
#
# Modules required: None
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_SPAdes.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_SCCmecFinder.sh sample_name	run_ID"
	echo "Output by default is sent to ${processed}/miseq_run_ID/sample_name/SCCmec"
	exit 0
elif [ -z "${2}" ]; then
	echo "Empty project id supplied to run_SPAdes.sh, exiting"
	exit 1
fi

if [[ ! -d "${processed}/${2}/${1}/SCCmec" ]]; then
	mkdir -p "${processed}/${2}/${1}/SCCmec"
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}/SCCmec"

python2 "${shareScript}/SCCmecFinder_v4.py"

#Script exited gracefully (unless something else inside failed)
exit 0
