#!/bin/bash -l

#$ -o run_sum.out
#$ -e run_sum.err
#$ -N run_sum
#$ -cwd
#$ -q short.q


#
# Description: Allows script to source this file to be able to pull out latest database filenames
#
# Usage ./get_latest_DBs.sh path_to_database_folder
#
# Output loction: std_out
#
# Modules required: None
#
# v1.0 (05/13/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./get_latest_DBs.sh path_to_database_folder"
	exit 0
elif [[ ! -d "${1}" ]]; then
	echo "Database location (${1}) does not exist, exiting"
	exit 2
fi

databases="${1}"
database_underscore_index=$(tr -dc '_' <<<"$databases" | awk '{ print length; }')
if [ -z ${database_underscore_index} ]; then
	database_underscore_index=0
fi
#echo $database_underscore_index
database_underscore_index=$(( database_underscore_index + 2 ))
#echo "Using ${1} as database location"

function get_ANI_REFSEQ {
	REFSEQ="NOT_FOUND"
	REFSEQ=$(find ${databases}/ANI -maxdepth 1 -name "REFSEQ_*.msh" -type f -printf '%p\n' | sort -k$database_underscore_index,$database_underscore_index -rt '_' -n | head -n 1)
	echo "${REFSEQ}"
}

function get_ANI_REFSEQ_Date {
	REFSEQ="NOT_FOUND"
	REFSEQ=$(find ${databases}/ANI -maxdepth 1 -name "REFSEQ_*.msh" -type f -printf '%p\n' | sort -k$database_underscore_index,$database_underscore_index -rt '_' -n | head -n 1)
	REFSEQ_date=$(echo "${REFSEQ}" | rev | cut -d'_' -f3 | rev )
	echo "${REFSEQ_date}"
}

function get_srst2 {
	ResGANNCBI_srst2="NOT_FOUND"
	ResGANNCBI_srst2=$(find ${databases}/star -maxdepth 1 -name "ResGANNCBI_*_srst2.fasta" -type f -printf '%p\n' | sort -k$database_underscore_index,$database_underscore_index -rt '_' -n | head -n 1)
	echo "${ResGANNCBI_srst2}"
}

function get_srst2_filename {
	ResGANNCBI_srst2="NOT_FOUND"
	ResGANNCBI_srst2=$(find ${databases}/star -maxdepth 1 -name "ResGANNCBI_*_srst2.fasta" -type f -printf '%p\n' | sort -k$database_underscore_index,$database_underscore_index -rt '_' -n | head -n 1)
	ResGANNCBI_srst2_filename=$(echo "${ResGANNCBI_srst2}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2)
	echo "${ResGANNCBI_srst2_filename}"
}

function get_ratio {
	NCBI_ratio="NOT_FOUND"
	NCBI_ratio=$(find ${databases}/ratio_DBs -maxdepth 1 -name "NCBI_Assembly_stats_*.txt" -type f -printf '%p\n' | sort -k4,4 -rt '_' -n | head -n 1)
	echo "${NCBI_ratio}"
}

function get_ratio_Date {
	NCBI_ratio="NOT_FOUND"
	NCBI_ratio=$(find ${databases}/ratio_DBs -maxdepth 1 -name "NCBI_Assembly_stats_*.txt" -type f -printf '%p\n' | sort -k4,4 -rt '_' -n | head -n 1)
	NCBI_ratio_date=$(echo "${NCBI_ratio}" | rev | cut -d'_' -f1 | cut -d'.' -f2 | rev)
	echo "${NCBI_ratio_date}"
}
