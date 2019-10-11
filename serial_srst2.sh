#!/bin/sh -l

#$ -o serial_srst2.out
#$ -e serial_srst2.err
#$ -N serial_srst2
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: A script to submit a list of isolates to the cluster to perform SRST2 AR on many isolates in parallel
#
# Usage: ./serial_srst2.sh path_to_list clobberness[keep|clobber]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/c-sstar/. Temp scripts will be in default_mass_qsubs_folder_from_config.sh/srst2_subs
#
# Modules required: None, run_srst2AR.sh will load srst2/0.2.0
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
	echo "Usage is ./serial_srst2.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions path_to_alt_database output_directory_for_scripts clobberness[keep|clobber]"
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

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to check all isolates on the list against the newest ResGANNCBI DB
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" = "clobber" ]]; then
		rm ${processed}/${project}/${sample}/srst2/${sample}__fullgenes__${ResGANNCBI_srst2_filename}_srst2__results.txt
		rm ${processed}/${project}/${sample}/srst2/${sample}__genes__${ResGANNCBI_srst2_filename}_srst2__results.txt
	fi
	echo ${counter}
	# Check if counter is below max number of concurrent submissions
	if [ ${counter} -lt ${max_subs} ]; then
		# Check if either one of the output files of srst2 files exist, skip submission if so
		if [[ ! -f "${processed}/${project}/${sample}/srst2/${sample}__genes__${ResGANNCBI_srst2_filename}_srst2__results.txt" ]] || [[ ! -f "${processed}/${project}/${sample}/srst2/${sample}__fullgenes__${ResGANNCBI_srst2_filename}_srst2__results.txt" ]]; then
			"${shareScript}/run_srst2AR.sh" "${sample}" "${project}"
		else
			echo "${project}/${sample} already has newest srst2 ResGANNCBI ${ResGANNCBI_srst2_filename}"
		fi
	fi
	counter=$(( counter + 1 ))
done

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "serial_srst2.sh has completed" "${global_end_time}" | mail -s "serial_srst2.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
