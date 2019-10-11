#!/bin/sh -l

#$ -o run_MLST.out
#$ -e run_MLST.err
#$ -N run_MLST
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
#  Description: Script to find mlst profile using tseeman's mlst script and pubmlst
# 	The most similar match is identified and provided for confirmation
#
# Usage: ./run_MLST.sh sample_name run_ID -f(force to certain DB) DB_to force_comparison_to
#
# Output location: default_config.sh_output_location/run_ID/sample_name/MLST/
#
# Modules required: mlst/2.16, perl/5.16.1-MT
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml mlst/2.16 perl/5.16.1-MT

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_MLST.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_MLST.sh   sample_name   run_ID [-f] [DB_to_force_to]"
	echo "Output is saved to ${processed}/run-id/sample_name/mlst"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id supplied to run_MLST.sh...exiting"
	exit 1
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}"

# Checks to see if a mlst folder already exists and creates it if not
if [ ! -d "$OUTDATADIR/MLST" ]; then
	echo "Creating $OUTDATADIR/MLST"
	mkdir -p "$OUTDATADIR/MLST"
fi

# Call mlst against the given mlst DB
#. "${shareScript}/module_changers/perl_5221_to_5161mt.sh"
if [[ "${3}" = "-f" ]]; then
	mlst --scheme "${4}" "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" > "${OUTDATADIR}/MLST/${1}_${4}.mlst"
# Call mlst using built autoidentifier
else
	mlst "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" > "${OUTDATADIR}/MLST/${1}.mlst"
fi
#. "${shareScript}/module_changers/perl_5161mt_to_5221.sh"
#Script exited gracefully (unless something else inside failed)

ml -mlst/2.16 -perl/5.16.1-MT

exit 0
