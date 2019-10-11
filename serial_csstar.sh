#!/bin/sh -l

#$ -o serial-cs.out
#$ -e serial-cs.err
#$ -N serial-cs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: A script to submit a list of isolates to the cluster to perform csstar on many isolates in parallel
#
# Usage: ./serial_csstar.sh path_to_list clobberness[keep|clobber] %ID(optional)[80|95|98|99|100]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/c-sstar/. Temp scripts will be in default_mass_qsubs_folder_from_config.sh/c-sstar_subs
#
# Modules required: None, run_c-sstar.sh will load Python3/3.5.2 and ncbi-blast+/LATEST
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
	echo "Usage is ./serial_csstar.sh path_to_list_file(single sample ID per line, e.g. B8VHY/1700128 (it must include project id also)) max_concurrent_submissions output_directory_for_scripts clobberness[keep|clobber] %ID(optional)[80|95|98|99|100]"
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

# Checks that value given for % Identity is one of the presets for csstar
if [[ "${3}" != 80 ]] && [[ "${3}" != 95 ]] && [[ "${3}" != 98 ]] && [[ "${3}" != 99 ]] && [[ "${3}" != 100 ]]; then
	echo "Identity is not one of the presets for csstar and therefore will fail, defaulting to 98..."
	sim="h"
	simnum=98
else
	if [ "${3}" == 98 ]; then
		sim="h"
	elif [ "${3}" == 80 ]; then
		sim="l"
	elif [ "${3}" == 99 ]; then
		sim="u"
	elif [ "${3}" == 95 ]; then
		sim="m"
	elif [ "${3}" == 100 ]; then
		sim="p"
	elif [ "${3}" == 40 ]; then
		sim="o"
	fi
	simnum=${3}
fi

# create an array of all samples in the list
arr=()
while IFS= read -r line || [ "$line" ];  do
	if [[ ! -z "${line}" ]]; then
		line=$(echo ${line} | tr -d '\n' | tr -d '\r')
		arr+=($line)
	fi
done < ${1}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"


# Create direcory to hold all temporary qsub scripts
counter=0

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit scripts to run default csstar on all samples on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" = "clobber" ]]; then
		echo "Removing old c-sstar for ${project}/${sample} and ${ResGANNCBI_srst2_filename}"
		rm ${processed}/${project}/${sample}/c-sstar/${sample}.${ResGANNCBI_srst2_filename}.gapped_${simnum}_sstar_summary.txt
		rm -r ${processed}/${project}/${sample}/c-sstar/${ResGANNCBI_srst2_filename}_gapped/
	fi
	#echo ${counter}-${project}-${sample}
	# Check if sample has a usable assembly file
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		#echo "Test"
		# Check if counter is below max number of concurrent submissions
		if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${ResGANNCBI_srst2_filename}.gapped_${simnum}_sstar_summary.txt" ]]; then
			"${shareScript}/run_c-sstar.sh" "${sample}" "g" "${sim}" "${project}"
		else
			echo "${project}/${sample} already has the newest ResGANNCBI (${ResGANNCBI_srst2_filename})"
		fi
	else
		echo "${project}/${sample} does not have an assembly to run csstar on"
	fi
	counter=$(( counter + 1 ))
done

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "serial_csstar.sh has completed" "${global_end_time}" | mail -s "serial_csstar.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
