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

run_tasks=11
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
run_name=$(basename ${run_to_check})
BAR_length=100
BAR_character='#'
BAR=$(printf "[%${BAR_length}s]" | tr ' ' $BAR_character)

# run_task AA
declare -A run_AA=( [1]="Copying Reads/Assemblies to project directory" [2]="Inverting list" [3]="Listing all isolates" [4]="Displaying isolates" [5]="Creating unique run Identifier" [6]="Catting list" [7]="running isolates" [8]="Creating bug array "[9]="Creating Seqlog" [10]="Creating run summary" [11]="Copying config and closing out run")
# isolate task AA
declare -A iso_AA=( [1]="Prepping FASTQ folder" [2]="Raw Read Quality count" [3]="BBDUK PhiX" [4]="Trimmomatic" [5]="Trimmed Read Quality Count" [6]="Kraken on reads" [7]="GOTTCHA" [8]="SRST2 AR" [9]="SPAdes Assembling" [10]="Trimming Assemmbly" [11]="Kraken on Assembly" [12]="16s Identification" [13]="Assembly QC" [14]="PROKKA" [15]="Rename Contig Headers" [16]="ANI" [17]="Taxon classification" [18]="BUSCO" [19]="c-SSTAR" [20]="GAMA" [21]="MLST" [22]="plasmidFinder" [23]="plasFlow" [24]="Check plasFlow assembly" [25]="c-SSTAR on plasFlow" [26]="plasmidFinder on PlasFlow" [27]="GAMA on plasFlow" [28]="Summarize isolate" [29]="Cleaning isolate")

#while 1; do
	pro_run_task_id=$(head -n1 ${run_to_check}/progress.txt | cut -d':' -f2)
	pro_Isolate_count=$(head -n2 ${run_to_check}/progress.txt | tail -n1 | cut -d':' -f2)
	current_Isolate_number=$(head -n3 ${run_to_check}/progress.txt | tail -n1 | cut -d':' -f2)
	isolate_index=$((current_Isolate_number + 1))
	current_Isolate_name=$(head -n${isolate_index} ${run_to_check}/${run_ID}_list.txt)
	pro_Isolate_task_number=$(tail -n1 ${run_to_check}/progress.txt | cut -d':' -f2)
	total_jobs=$(( run_tasks + pro_Isolate_count * tasks_per_isolate ))
	current_Isolate_progress=$(( 100 * pro_Isolate_task_number / tasks_per_isolate ))
	jobs_completed=$(( current_Isolate_number * tasks_per_isolate - tasks_per_isolate + pro_run_task_id + pro_Isolate_task_number))
	total_progress=$(( 100 * jobs_completed / total_jobs ))
	echo -e "${pro_run_task_id}	${pro_Isolate_count}	${current_Isolate_number}	${pro_Isolate_task_number}	${total_jobs}	${jobs_completed}\n\n\n"
	echo "${current_Isolate_progress}"
	echo "${total_progress}"
	echo -ne "\r${BAR:0:$current_Isolate_progress}(${current_Isolate_progress}%-${current_Isolate_name}-${iso_AA[${pro_Isolate_task_number}]})"
	echo -ne "\r${BAR:0:$total_progress}(${total_progress}%-${run_AA[${pro_run_task_id}]})"
	sleep 1
	if [[ "${total_progress}" -eq 100 ]]; then
		echo "Run is complete!!!"
		exit
	fi
#done
