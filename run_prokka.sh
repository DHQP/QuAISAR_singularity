#!/bin/sh -l

#$ -o run_prokka.out
#$ -e run_prokka.err
#$ -N run_prokka
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh
#. "${mod_changers}/perl_5221_to_5123.sh"

#
# Description: Runs prokka gene identifier on sample to discover all identifiable genes. Also necessary for downstream busco processing
#
# Usage ./run_prokka.sh   sample_name   run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/prokka
#
# Modules required: prokka/1.12, perl/5.12.3
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml prokka/1.12 perl/5.12.3

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_prokka.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_prokka.sh sample_name run_ID"
	echo "Output is saved to ${processed}/miseq_run_ID/sample_name/prokka"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id supplied to run_prokka.sh, exiting"
	exit 1
fi

# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${2}/${1}"

# Checks for existence of prokka output folder and deletes and recreates it if there
if [ -d "$OUTDATADIR/prokka" ]; then  #removes old prokka results before continuing (it will complain otherwise)
	echo "Removing old prokka results $OUTDATADIR/prokka"
	rm -rf "$OUTDATADIR/prokka"
fi

### Prokka to identify genes in ###
echo "Running Prokka for gene identification"
# Run prokka
prokka --outdir "${OUTDATADIR}/prokka" "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta"
#Rename all PROKKA files with sample name instead of date
if [ ! -d "${OUTDATADIR}/prokka" ]; then
	echo "prokka did not complete...exiting"
	# reload perl to 5.22.1 before exiting
	. "${mod_changers}/perl_5123_to_5221.sh"
	exit 1
fi
#echo "About to rename files"
for pfile in ${OUTDATADIR}/prokka/*.*; do
	fullname=$(basename "${pfile}")
	ext="${fullname##*.}"
	echo "Renaming ${pfile} to ${OUTDATADIR}/prokka/${1}_PROKKA.${ext}"
	mv "${pfile}" "${OUTDATADIR}/prokka/${1}_PROKKA.${ext}"
done

#Script exited gracefully (unless something else inside failed)

ml -prokka/1.12 -perl/5.12.3

exit 0
