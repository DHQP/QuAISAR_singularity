#!/bin/sh -l

#$ -o serial_GAMA.out
#$ -e serial_GAMA.err
#$ -N serial_GAMA
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: A script to submit a list of isolates to the cluster to perform GAMA on many isolates in parallel
#
# Usage: ./serial_GAMA.sh path_to_list clobberness[keep|clobber]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/GAMA/. Temp scripts will be in default_mass_qsubs_folder_from_config.sh/GAMA_subs
#
# Modules required: None, run_GAMA.sh will load Python3/3.5.2 and blat/35
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./serial_GAMA.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions path_to_alt_database output_directory_for_scripts clobberness[keep|clobber]"
	exit 1
elif [[ ! -f "${1}" ]]; then
	echo "${1} (list) does not exist...exiting"
	exit 1
elif [[ -z "${2}" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${2}" != "keep" ]] && [[ "${2}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 5th parameter...exiting"
	exit 1
else
	clobberness="${2}"
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [ "$line" ];  do
  arr+=("$line")
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

# Create counter and set max number of concurrent submissions
counter=0

"${shareScript}/clean_list.sh" "${1}"

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to check all isolates on the list against the newest ResGANNCBI DB
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" = "clobber" ]]; then
		rm ${processed}/${project}/${sample}/GAMA/${sample}_${ResGANNCBI_srst2_filename}.GAMA
	fi
	echo ${counter}
	# Check if counter is below max number of concurrent submissions
	if [ ${counter} -lt ${max_subs} ]; then
		# Check if the output file of GAMA exist, skip submission if so
		if [[ ! -f "${processed}/${project}/${sample}/GAMA/${sample}_${ResGANNCBI_srst2_filename}.GAMA" ]]; then
			"${shareScript}/run_GAMA.sh" "${sample}" "${project}" "-c"
		else
			echo "${project}/${sample} already has newest GAMA ResGANNCBI ${ResGANNCBI_srst2_filename}"
		fi
	fi
	counter=$(( counter + 1 ))
done


echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "serial_GAMA.sh has completed" "${global_end_time}" | mail -s "serial_GAMA.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
