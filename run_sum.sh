#!/bin/sh -l

#$ -o run_sum.out
#$ -e run_sum.err
#$ -N run_sum
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Creates a summary file for the run and prints out a one word status and a short description of each step being reported
#
# Usage ./run_sum.sh run_ID
#
# Output loction: default_config.sh_output_location/run_ID/
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
	echo "Empty run_ID supplied to run_sum.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_sum.sh miseq_run_ID -vo(optional)"
	echo "Output is saved to ${processed}/miseq_run_ID"
	exit 0
fi

echo "Checking for ${processed}/${1}/${1}_list(_ordered).txt"

# Checks for existence of list files in specific order
if [[ -z ${2} ]]; then
	if [[ -f ${processed}/${1}/${1}_list_ordered.txt ]]; then
		list="${processed}/${1}/${1}_list_ordered.txt"
	elif [[ -f ${processed}/${1}/${1}_list.txt ]]; then
		list="${processed}/${1}/${1}_list.txt"
	else
		echo "No list file exists, cannot do a summary, unless I add in an automagic creator later"
		exit
	fi
	type="project"
else
	type="list"
	list=${1}
fi

# Gets todays date to show when summary was run
runsumdate=$(date "+%Y_%m_%d_at_%Hh_%Mm")
echo "Creating run summary at ${runsumdate}"
# Status of each individual sample is updated in its own folder and the run_summary file
if [[ "${type}" = "project" ]]; then
	sum_name="${1}_run_summary_at_${runsumdate}.sum"
	echo "named as project"
else
	sum_name="list_summary_at_${runsumdate}.sum"
	echo "named as list"
fi

# Run validate_piperun.sh on every sample in the list and cat output into one summary file
while IFS= read -r samples || [ -n "$samples" ]; do
	echo ${file}
	file=$(echo "${samples}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	proj=$(echo "${samples}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	"${shareScript}/validate_piperun.sh" "${file}" "${proj}" > "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"
	if [[ "${type}" = "project" ]]; then
		cat "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" >> "${processed}/${proj}/${sum_name}"
	else
		cat "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" >> "${3}/${sum_name}"
	fi
done < ${list}
