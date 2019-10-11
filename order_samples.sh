#!/bin/sh -l

#$ -o order_samples.out
#$ -e order_samples.err
#$ -N order_samples
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh

#
# Description: Creates a txt list file that contains run and samples for a given MiSeq run that matches the order of the output for the MMB_Seq log
#
# Usage: ./order_samples.sh MiSeq_Run_ID
#
# Output location: /deafult_config.sh_output_location/Run_ID/Run_ID_list_ordered.txt
#
# Modules required: Python2/2.7.13
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml Python2/2.7.13

#  Function to print out help blurb
show_help () {
	echo "Usage is ./order_samples.sh -p project_name"
	echo "Output is saved to path_to_folder"
}

# Parse through command line options
options_found=0
while getopts ":h?n:p:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			project=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

#### #copy MMBSEQLog to local
> "${processed}/${project}/${project}_list_ordered.txt"

echo "${processed}/${project}/2019_MMBSeq_Log.xlsx"

# Copy the newest log file to the local directory
cp "${local_DBs}/Seqlog_copies/2019_MMBSeq_Log.xlsx" "${processed}/${project}/2019_MMBSeq_Log.xlsx"

# Convert log file to csv format for searchability
python3 ${shareScript}/xlsx_converter.py "${processed}/${project}/2019_MMBSeq_Log.xlsx" "FY19 Miseq Isolate Log" > "${processed}/${project}/2019_MMBSeq_Log.tsv"

echo "Excel file: 2019_MMBSeq_Log.xlsx has been converted to TSV"

# Parse log file csv until run_if matches
counter=0
while IFS= read -r var || [ -n "$var" ]; do
	# Check the format of the city/state column in the log file to determine how many tabs need to be used to find run_ID in line
	#echo "checking ${var}"
	# city_state=$(echo "${var}" | cut -d',' -f17)
	# echo "|^| "${#city_state}" : "${city_state}
	# if [[ ${#city_state} -eq 2 ]] || [[ ${city_state} = "unknown" ]] || [[ ${city_state} = "Unknown" ]] || [[ ${city_state} = "UNKNOWN" ]] || [[ "${city_state}" = "" ]]; then
	# 	echo "Doing if"
	# 	line_project=$(echo "${var}" | cut -d',' -f21)
	# else
	# 	echo "doing else"
	# 	line_project=$(echo "${var}" | cut -d',' -f22)
	# fi
	line_project=$(echo "${var}" | cut -d'	' -f21)
	# echo "${line_project}:${project}"
	# If the run_ID matches, then add ID to list (automatically placing them in the proper order)
	if [[ "${line_project}" = "${project}" ]]; then
		line_id=$(echo "${var}" | cut -d'	' -f3)
		#echo "Adding ${counter}: ${project}/${line_id}"
		echo "${project}/${line_id}" >> "${processed}/${project}/${project}_list_ordered.txt"
	else
		#echo "${counter}-${line_project} not in ${project}"
		:
	fi
	counter=$(( counter + 1 ))
done < ${processed}/${project}/2019_MMBSeq_Log.tsv


# Remove intermediate files from sorting
rm -r ${processed}/${project}/2019_MMBSeq_Log.tsv
rm -r ${processed}/${project}/2019_MMBSeq_Log.xlsx

# Check if the sorted file has content, else delete it since something went wrong
if [[ ! -s "${processed}/${project}/${project}_list_ordered.txt" ]]; then
	echo "Isolates were not able to be sorted, something wrong with MiSeq Log entries or list file, or....?"
	rm -r ${processed}/${project}/${project}_list_ordered.txt
	exit
else
	echo "sorted file contains entries"
	:
fi

ml -Python2/2.7.13
