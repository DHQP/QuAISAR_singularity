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
# Description: Helper for the main quaisar script to allow easy visualization of progress of run
#
# Usage ./quaisar-progress.sh full_path_to_run_ID
#
# Output loction: screen
#
# Modules required: None
#
# v1.0 (03/19/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

run_tasks=10
tasks_per_isolate=29

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty run_ID supplied to quaisar-progress.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_sum.sh path_to_run_folder"
	echo "Output is saved to ${output_dir}/miseq_run_ID"
	exit 0
fi

run_to_check=${1}
while 1; do
	pro_run_task_id=$(head -n1 ${run_to_check}/progress.txt | cut -d':' -f2)
	pro_Isolate_count=$(head -n2 ${run_to_check}/progress.txt | tail -n1 | cut -d':' -f2)
	current_Isolate_number=$(head -n3 ${run_to_check}/progress.txt | tail -n1 | cut -d':' -f2)
	pro_Isolate_task_number=$(tail -n1 ${run_to_check}/progress.txt | cut -d':' -f2)
	total_jobs=$(( run_tasks + pro_Isolate_count * tasks_per_isolate ))
	current_Isolate_progress=$(( 100 * pro_Isolate_task_number / tasks_per_isolate ))
	jobs_completed=$(( current_Isolate_number * tasks_per_isolate - tasks_per_isolate + pro_run_task_id + pro_Isolate_task_number))
	total_progress=$(( 100 * jobs_completed / total_jobs ))
	echo -e "${pro_run_task_id}	${pro_Isolate_count}	${current_Isolate_number}	${pro_Isolate_task_number}	${total_jobs}	${jobs_completed}\n\n\n"
	echo "${current_Isolate_progress}"
	echo "${total_progress}"
	sleep 1
	if [[ "${total_progress}" -eq 100 ]];
		echo "Run is complete!!!"
		exit
	fi
done
