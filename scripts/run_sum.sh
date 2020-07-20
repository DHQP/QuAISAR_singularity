#!/bin/sh -l

#$ -o run_sum.out
#$ -e run_sum.err
#$ -N run_sum
#$ -cwd
#$ -q short.q

#
# Description: Creates a summary file for the run and prints out a one word status and a short description of each step being reported
#
# Usage ./run_sum.sh path_to_run_folder path_to_scripts_folder path_to_database_folder
#
#
# Output loction: path_to_run_folder/run_name_run_summary_at_time.sum
#
# Modules required: None
#
# v1.0b (05/14/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_sum.sh path_to_run_folder path_to_scripts_folder path_to_database_folder"
	echo "Output is saved to path_to_run_folder/run_name_run_summary_at_time.sum"
	exit 0
elif [[ ! -d "${1}" ]]; then
	echo "Path ($1) does not exit, exiting run_sum.sh"
	exit 2
elif [[ -z "${2}" ]]; then
	echo "Empty script path supplied, exiting run_sum.sh"
	exit 3
elif [[ ! -d "${2}" ]]; then
	echo "Script path ($2) does not exsit, exiting run_sum.sh"
	exit 4
elif [[ ! -f "${2}/validate_piperun.sh" ]]; then
	echo "validate_piperun.sh does not exist in $2, exiting run_sum.sh"
	exit 5
elif [[ -z "${3}" ]]; then
	echo "Empty database pathsupplied, exiting run_sum.sh"
	exit 6
elif [[ ! -d "${3}" ]]; then
	echo "Database path ($3) does not exsit, exiting run_sum.sh"
	exit 7
fi

OUTDATADIR="${1}"
# Based upon standard naming protocols pulling 2nd to last portion of path off should result in proper project name
project_name=$(echo "${OUTDATADIR}" | rev | cut -d'/' -f1 | rev)
databases=${2}
scripts="${3}"

echo "Checking for ${OUTDATADIR}/${project_name}/${project_name}_list(_ordered).txt"
list=$(find ${OUTDATADIR} -name "*_list_ordered.txt")
if [[ ! -f "${list}" ]]; then
	list=$(find ${OUTDATADIR} -name "*_list.txt")
fi
if [[ ! -f "${list}" ]]; then
	echo "No list found (${OUTDATADIR}/*_list(_ordered).txt), exiting"
	exit
fi

# Gets todays date to show when summary was run
runsumdate=$(date "+%Y_%m_%d_at_%Hh_%Mm")
echo "Creating run summary at ${runsumdate}"

# Run validate_piperun.sh on every sample in the list and cat output into one summary file
while IFS= read -r samples || [ -n "$samples" ]; do
	sample_name=$(echo "${samples}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project_name_internal=$(echo "${samples}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	if [[ "${project_name}" == "${project_name_internal}" ]]; then
		"${3}/validate_piperun.sh" "${OUTDATADIR}/${sample_name}" "${databases}" "${scripts}" > "${OUTDATADIR}/${sample_name}/${sample_name}_pipeline_stats.txt"
	fi
done < ${list}
