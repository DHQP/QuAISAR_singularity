${local_DBs}/singularities/srst2.simg#!/bin/bash -l

#$ -o quaisar_X.out
#$ -e quaisar_X.err
#$ -N quasX
#$ -cwd
#$ -q all.q

#
# Description: The full QuAISAR-H pipeline start to end serially with singularity containers
#
# Usage: ./quaisar_singularity.sh -i location_of_reads 1|2|3|4 -o name_of_output_folder -p project_name [-s full_path_to_script_folder] [-r] [-d full_path_to_database_folder] [-c config.sh full_path_to_config_file]"
#		filename postfix numbers are as follows 1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz 5: Asssemblies (.fasta)"
#		Reads can be gzipped or raw, but if your files are not named in any one of these formats, they will need to be renamed before running them through the pipeline
#		If you are submitting assemblies, use 1 as the value
#
# Output location: A folder with the name given for the -p flag will be created under the folder given with the -o flag (/output/project_name)
#
# v1.1 (11/17/2021)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Initialize progress variables
run_location=""
run_task_id=0
task_number=0
isolate_number=0
isolate_count=0

version_type="Quaisar-Singularity"
version_num="qs1.1"

#ml singularity
singularity --version

# Will be called throughout the script to write current progress for inquisitive minds and to restart run from where it was murdered
# parameters need to be as follows
# 1 - path_to_run_folder
# 2 - run_task_id
# 3 - total isolate count
# 4 - isolate number (out of however many on the list)
# 5 - isolate task number
function write_Progress() {
	echo -e "run_task_ID:${run_task_id}\nTotal_isolates:${isolate_count}\nCompleted_isolates:${isolate_number}\nCurrent_isolate_task_number:${task_number}" > ${PROJDATADIR}/progress.txt
}

# Checking for proper number of arguments from command line
if [[ $# -lt 1  || $# -gt 13 ]]; then
	echo -e "\\n\\n\\n"
	echo -e "Usage: ./quaisar_singularity.sh -i location_of_reads -p project_name [-o name_of_output_folder] [-s full_path_to_script_folder] [-r] [-a] [-d full_path_to_database_folder] [-c config.sh full_path_to_config_file]"
	echo -e "Reads filenames need to have a postfix in one of the following _S*_L001_R*_00*.fastq[.gz], _S*_R*_0*X.fastq[.gz], _RX_00*.fastq[.gz], _[R]*.fastq[.gz]."
	echo -e "Assembly filenames need to have a postfix of .fasta or .fna"
	echo -e "If your reads are not named in any one of these formats, they will need to be renamed before running them through the pipeline"
	echo -e "Reads can be gzipped or raw, but if you are submitting assemblies, use 1 as the value"
	echo -e "Additional functions/flags: \n\t -s If you would like to reference and run pipeline scripts installed in an alternate location \n\t -r if you would like to retry the list of samples if they failed during assembly \n\t -d if you would like to reference a different location for databases, but must contain all necessary for pipeline \n\t -a Run pipeline from assemblies"
	echo -e "\\n\\n\\n"
fi

# Checks the arguments and sets some default variables
nopts=$#
do_download="false"
global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
requestor=$(whoami)
PROJECT="${requestor}_${global_time}"
assemblies="false"
current_directory=$(pwd)
config_file=${current_directory}/config.sh
. ${config_file}

for ((i=1 ; i <= nopts ; i++)); do
	#echo "${1} ${2}"
	case "${1}" in
		#Help/Usage section
		-h | --help)
			echo -e "\\n\\n\\n"
			echo -e "Usage: ./quaisar_singularity.sh -i location_of_reads -o name_of_output_folder -p project_name [-s full_path_to_script_folder] [-r] [-a] [-d full_path_to_database_folder] [-c config.sh full_path_to_config_file]"
			echo -e "Reads filenames need to have a postfix in one of the following _S*_L001_R*_00*.fastq[.gz], _S*_R*_0*X.fastq[.gz], _RX_00*.fastq[.gz], _[R]*.fastq[.gz]."
			echo -e "Assembly filenames need to have a postfix of .fasta or .fna"
			echo -e "If your reads are not named in any one of these formats, they will need to be renamed before running them through the pipeline"
			echo -e "Reads can be gzipped or raw, but if you are submitting assemblies, use 1 as the value"
			echo -e "Additional functions/flags: \n\t -s If you would like to reference and run pipeline scripts installed in an alternate location \n\t -r if you would like to retry the list of samples if they failed during assembly \n\t -d if you would like to reference a different location for databases, but must contain all necessary for pipeline \n\t -a Run pipeline from assemblies"
			echo -e "\\n\\n\\n"
			exit 0
			;;
		#Import the config file to set other minor locations to be used within the script
		-c | --config)
			config_file="$2"
			if [ -f "${config_file}" ]; then
				. "${config_file}"
				BASEDIR="${output_dir}"
			else
				echo "Can not find config file, $2"
				exit 22
			fi
			shift 2
			;;
		#Gets name of folder that FASTA files will be in
		-i | --in-dir)
			INDATADIR="$2"
			if [[ -d  ${INDATADIR} ]]; then
				do_download="true"
				if [ "${INDATADIR:0:1}" != "/" ] && [ "${INDATADIR:0:1}" != "$" ]; then
					echo "${INDATADIR}"
					echo 'ERROR: The full input path was not specified.' >&2
					exit 1
				fi
			else
					echo "FASTQ folder ${INDATADIR} does not exist...exiting"
					exit 1
			fi
			indir_set="true"
			#is_full_run="false"
			#echo "$INDATADIR $2"
			shift 2
			;;
		#Gets output directory name of folder that all output files will be stored
		-o | --out-dir)
			echo "Setting output directory as: ${2}"
			output_dir="$2"
			shift 2
			;;
			#Gets output directory name of folder that all output files will be stored
		-p | --project_name)
			echo "Setting project name as: ${2}"
			PROJECT="$2"
			shift 2
			;;
		-d | --database)
			echo "Setting database location as: ${2}"
			local_DBs="$2"
			shift 2
			;;
		-s | --scripts_location)
			echo "Setting script location as: ${2}"
			src="$2"
			shift 2
			;;
		-r | --retry_from_assembly)
			assemblies="retry"
			shift
			;;
		-a | --assemblies)
			assemblies="true"
			;;
		#Captures any other characters in the args
		\?)
			echo "ERROR: ${BOLD}$2${NORM} is not a valid argument" >&2
			usage
			exit 1
			;;
	esac
done


# Check if specific conda environment exists
# Havent found command yet

# Turn on Conda environment?
# figure out best practices install location
# conda env create -f ~/py36_biopython/py36_biopython_singularity.yml

prereqs="true"
missing_names=()

conda activate py36_biopython

echo "Checking for dependencies and databases"
# Check for required software (python3 and singularity)
conda_call_lines=$(conda list | wc -l)
singularity_version=$(singularity --version | cut -d' ' -f3 | cut -d'.' -f1)
singularity_release=$(singularity --version | cut -d' ' -f3)
python_release=$(python -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
python_version=$(echo "${python_release}" | cut -d'.' -f1)

if [[ "${conda_call_lines}" -gt 1 ]]; then
	:
else
	missing_names=("${missing_names[@]}" conda)
fi

if [[ "${python_version}" = "3" ]]; then
	python_command="python"
	#echo "Python $python_release is installed, please continue"
elif [[ "${python_version}" = "2" ]]; then
	python_version=$(python3 --version | cut -d' ' -f2 | cut -d'.' -f1)
	if [[ "${python_version}" = "3" ]]; then
		python_command="python3"
		python_release=$(python3 --version | cut -d' ' -f2)
		#echo "Python $python_release is installed, please continue"
	else
		#echo -e "\nPython3.x not installed, can not proceed\n"
		prereqs="false"
		missing_names=("${missing_names[@]}" python3)
	fi
else
	echo -e "\nPython3.x not installed, can not proceed\n"
	echo -e "Python version =${python_release}"
	prereqs="false"
	missing_names=("${missing_names[@]}" python3)
fi


echo "input - ${INDATADIR}"
echo "output - ${output_dir}"
echo "project - ${PROJECT}"
echo "dbs - ${local_DBs}"
echo "src - ${src}"


#echo "${singularity_version}-${singularity_release}:${python_version}-${python_release}:${python_command}"


# Check that biopython is installed
${python_command} -c "import Bio"
bio_installed=$(echo $?)
if [[ "${bio_installed}" -eq 0 ]]; then
	#echo "Biopython is installed, please continue"
	:
else
	#echo -e "\nBiopython not installed, can not proceed\n"
	prereqs="false"
	missing_names=("${missing_names[@]}" biopython)
fi

# Check singularity version...dirty fix for dirty until i can figure out why it doesnt give correct version #
if [[ "${singularity_version}" -ge 3 ]] || [[ "${singularity_version}" = "d5eaf8a+dirty" ]]; then
	#echo "Singularity ${singularity_release} is installed, please continue"
	:
else
	#echo -e "\nSingularity 3.x(+) is not installed, can not continue\n"
	prereqs="false"
	missing_names=("${missing_names[@]}" singularity)
fi

db_output=$(${src}/database_checker.sh ${local_DBs} | tail -n1)
num_missing=$(echo "${db_output}" | cut -d' ' -f3)
db_missing=$(echo "${db_output}" | cut -d'(' -f2 | cut -d')' -f1)

if [[ "${num_missing}" -gt 0 ]]; then
	prereqs="false"
	missing_names=("${missing_names[@]}" "${db_missing[@]}")
fi

if [[ "${prereqs}" = "false" ]]; then
	echo "A dependency or database is missing:"
	for missing in "${missing_names[@]}"; do
		echo -e "${missing}\n"
	done
	exit
else
	echo "Everything is preinstalled and ready to run"
fi


list_path="${output_dir}/${PROJECT}/${PROJECT}_list.txt"
if [[ ! -d ${output_dir} ]]; then
	mkdir -p ${output_dir}
fi

# Set database names to use
. "${src}/get_latest_DBs.sh" "${local_DBs}"
ResGANNCBI_srst2=$(get_srst2)
ResGANNCBI_srst2_filename=$(get_srst2_filename)
REFSEQ=$(get_ANI_REFSEQ)
REFSEQ_date=$(get_ANI_REFSEQ_Date)
NCBI_ratio=$(get_ratio)
NCBI_ratio_date=$(get_ratio_Date)

if [[ -z "${ResGANNCBI_srst2_filename}" ]]; then
	echo "ResGANNCBI_srst2_filename (db version) is empty, cant perform c-SSTAR, srst2, or GAMMA AR analysis"
fi

if [[ -z "${ResGANNCBI_srst2}" ]]; then
	echo "ResGANNCBI_srst2_filename is empty, cant perform c-SSTAR, srst2, or GAMMA AR analysis"
fi

if [[ -z "${REFSEQ_date}" ]]; then
	echo "REFSEQ_date (db version) is empty, cant perform c-SSTAR, srst2, or GAMMA AR analysis"
fi

if [[ -z "${REFSEQ}" ]]; then
	echo "REFSEQ ANI Database is empty, cant perform ANI"
fi


# Short print out summary of run settings
echo -e "Source folder: ${INDATADIR}\\nOutput folder: ${output_dir}\\nList based analysis:  ${list_path}"

# Sets folder to where files will be downloaded to
PROJDATADIR="${output_dir}/${PROJECT}"
if [ ! -d "${PROJDATADIR}" ]; then
	echo "Creating $PROJDATADIR"
	mkdir -p "${PROJDATADIR}"
fi
if [ -f "${PROJDATADIR}/${PROJECT}_list.txt" ]; then
	mv "${PROJDATADIR}/${PROJECT}_list.txt" "${PROJDATADIR}/${PROJECT}_list_original.txt"
fi

# Task: Copies reads/assemblies from source location to working directory and creates a list of IDs
write_Progress
run_task_id=1
if [[ "${assemblies}" == "true" ]]; then
	# Goes through given Assemblies folder
	echo "${INDATADIR}"
	for file in ${INDATADIR}/*
	do
		# Check if file is a recognized assembly format extension
		if [[ "${file}" = *.fasta ]] || [[ "${file}" = *.fna ]]; then
			isolate_name=$(basename -- "$file")
			extension="${isolate_name##*.}"
			sample="${isolate_name%.*}"

			mkdir -p ${PROJDATADIR}/${sample}/Assembly
			cp ${file} ${PROJDATADIR}/${sample}/Assembly/scaffolds.fasta
			echo -e "${PROJECT}/${sample}" >> "${PROJDATADIR}/${PROJECT}_list.txt"
		else
			echo "${file} is not an fna or fasta file, not acting on it"
		fi
		# Invert list so that the important isolates (for us at least) get run first
		if [[ -f "${PROJDATADIR}/${PROJECT}_list.txt" ]]; then
			sort -k2,2 -t'/' -r "${PROJDATADIR}/${PROJECT}_list.txt" -o "${PROJDATADIR}/${PROJECT}_list.txt"
		fi
	done
else
	# Goes through given reads folder
	echo "${INDATADIR}"
	for file in ${INDATADIR}/*
	do
		# Check if file is a zipped reads file
		if [[ "${file}" = *.fastq.gz ]] || [[ "${file}" = *.fastq ]] && [[ "${file}" != *_L001_I1_001.fastq.gz ]]; then
			full_sample_name=${file##*/}
			echo ${full_sample_name}
			# Extracts filename keeping only isolate ID, if it matches standard miseq naming
			if [[ ${full_sample_name} =~ _S[0-9]+_L[0-9]+_R[1|2]_00[0-9]+\.fast.+$ ]]; then #*"_S"*"_L001_RX_00X.fastq.gz" ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f5- | rev)
				if [[ "${full_sample_name}" = *"_R1_"* ]]; then
					full_sample_name_pair=${full_sample_name/_R1_/_R2_}
					current_read=1
				elif [[ "${full_sample_name}" = *"_R2_"* ]]; then
					full_sample_name_pair="${full_sample_name/_R2_/_R1_}"
					current_read=2
				fi
				# postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2,3 | rev)

			elif [[ ${full_sample_name} =~ _S[0-9]+_R[1|2]_00[0-9]+\.fast.+$ ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f4- | rev)
				if [[ "${full_sample_name}" = *"_R1_"* ]]; then
					full_sample_name_pair=${full_sample_name/_R1_/_R2_}
					current_read=1
				elif [[ "${full_sample_name}" = *"_R2_"* ]]; then
					full_sample_name_pair="${full_sample_name/_R2_/_R1_}"
					current_read=2
				fi
				#postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2,3 | rev)

			elif [[ ${full_sample_name} =~ _R[1|2]_00[0-9]+\.fast.+$ ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f3- | rev)
				if [[ "${full_sample_name}" = *"_R1_"* ]]; then
					full_sample_name_pair=${full_sample_name/_R1_/_R2_}
					current_read=1
				elif [[ "${full_sample_name}" = *"_R2_"* ]]; then
					full_sample_name_pair="${full_sample_name/_R2_/_R1_}"
					current_read=2
				fi
				#postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2 | rev)

			elif [[ ${full_sample_name} =~ _R[1|2]\.fast.+$ ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f2- | rev)
				if [[ "${full_sample_name}" = *"_R1.fast"* ]]; then
					full_sample_name_pair=${full_sample_name/_R1.fast/_R2.fast}
					current_read=1
				elif [[ "${full_sample_name}" = *"_R2.fast"* ]]; then
					full_sample_name_pair="${full_sample_name/_R2.fast/_R1.fast}"
					current_read=2
				fi
				#postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1 | rev)

			elif [[ ${full_sample_name} =~ _[1|2]\.fast.+$ ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f2- | rev)
				if [[ "${full_sample_name}" = *"_1.fast"* ]]; then
					full_sample_name_pair=${full_sample_name/_1.fast/_2.fast}
					current_read=1
				elif [[ "${full_sample_name}" = *"_2.fast"* ]]; then
					full_sample_name_pair="${full_sample_name/_2.fast/_1.fast}"
					current_read=2
				fi
				#postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1 | rev)

			else
				echo "Does not match any of the expected filenameing conventions - ${full_sample_name}"
				exit
			fi

			#long_name=$(echo "${full_sample_name}" | cut -d'_' -f1,2,3)
			echo "Short: ${short_name}"
			#echo "Does ${full_sample_name} match *${match}"
			source_path=$(dirname "${file}")

			# Skip file if it happens to be undetermined
	    if [[ "${short_name}" == "Undetermined" ]]; then
				echo "found undetermined (${file})"
				continue
			# If the file matches the postfix given in the arguments proceed with moving and unzipping to the output directory
			else
				# Creates output folder
				if [ ! -d "${PROJDATADIR}/${short_name}" ]; then
					echo "Creating $PROJDATADIR/${short_name}"
					mkdir -p "${PROJDATADIR}/${short_name}"
					echo "Creating ${PROJDATADIR}/${short_name}/FASTQs"
					mkdir -p "${PROJDATADIR}/${short_name}/FASTQs"
				fi
				# Announces name of file being unzipped and then unzips it to the FASTQs folder for the matching sample name. Files are shortened to just name_R1_001.fastq or name_R2_001.fastq
				echo "Retrieving ${source_path}/${full_sample_name} and ${full_sample_name_pair}"
				#if [[ "${match}" -eq 1 ]] || [[ "${match}" -eq 4 ]]; then
				if [[ "${current_read}" -eq 1 ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						if [[ -f "${source_path}/${full_sample_name_pair}" ]]; then
							if [[ -f "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" ]] && [[ -f "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" ]]; then
								echo "${short_name} already has both zipped FASTQs (Probably done when R2 was found, this is the R1 tester)"
							else
								echo "Moving ${source_path}/${full_sample_name} and ${source_path}/${full_sample_name_pair} to ${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz and ${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								cp "${source_path}/${full_sample_name}" "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
								cp "${source_path}/${full_sample_name_pair}" "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
							fi
						else
							echo "No R2 found for ${source_path}/${full_sample_name}"
							cp "${source_path}/${full_sample_name}" "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						fi
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						if [[ -f "${source_path}/${full_sample_name_pair}" ]]; then
							if [[ -f "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" ]] && [[ -f "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" ]]; then
								echo "${short_name} already has both unzipped FASTQs (Probably done when R2 was found, this is the R1 tester)"
							else
								echo "Zipping and not clumping ${source_path}/${full_sample_name} and ${source_path}/${full_sample_name_pair} to ${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz and ${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								gzip -c "${source_path}/${full_sample_name}" > "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
								gzip -c "${source_path}/${full_sample_name_pair}" > "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
							fi
						else
							echo "No R2 found for ${source_path}/${full_sample_name}"
							gzip -c "${source_path}/${full_sample_name}" > "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						fi
					fi
				elif [[ "${current_read}" -eq 2 ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						if [[ -f "${source_path}/${full_sample_name_pair}" ]]; then
							if [[ -f "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" ]] && [[ -f "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" ]]; then
								echo "${short_name} already has both zipped FASTQs (Probably done when R1 was found, this is the R2 tester)"
							else
								echo "Not Clumping ${source_path}/${full_sample_name} and ${source_path}/${full_sample_name_pair} to ${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz and ${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								cp "${source_path}/${full_sample_name}" "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								cp "${source_path}/${full_sample_name_pair}" "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
							fi
						else
							echo "No R1 found for ${source_path}/${full_sample_name}"
							cp "${source_path}/${full_sample_name}" "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
						fi
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						if [[ -f "${source_path}/${full_sample_name_pair}" ]]; then
							if [[ -f "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" ]] && [[ -f "${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" ]]; then
								echo "${short_name} already has both zipped FASTQs (Probably done when R1 was found, this is the R2 tester)"
							else
								echo "Zipping and not clumping ${source_path}/${full_sample_name} and ${source_path}/${full_sample_name_pair} to ${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz and ${PROJDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								gzip -c "${source_path}/${full_sample_name}" > "${PROJDATADIR}/${short_name}/FASTQs/temp/${short_name}_R1_001.fastq.gz"
								gzip -c "${source_path}/${full_sample_name_pair}" > "${PROJDATADIR}/${short_name}/FASTQs/temp/${short_name}_R2_001.fastq.gz"
							fi
						else
							echo "No R1 found for ${source_path}/${full_sample_name}"
							gzip -c "${source_path}/${full_sample_name}" > "${PROJDATADIR}/${short_name}/FASTQs/temp/${short_name}_R2_001.fastq.gz"
						fi
					else
						"Current read = ${current_read} (hint: its not 1 or 2), so it should never get here anyway"
					fi
					if [[ -f "${PROJDATADIR}/${PROJECT}_list.txt" ]]; then
						if grep -Fxq "${PROJECT}/${short_name}" "${PROJDATADIR}/${PROJECT}_list.txt"
						then
							echo -e "${PROJECT}/${short_name} already on list ${PROJDATADIR}/${PROJECT}_list.txt, not adding again"
						else
							echo -e "${PROJECT}/${short_name}" >> "${PROJDATADIR}/${PROJECT}_list.txt"
						fi
					else
						echo -e "${PROJECT}/${short_name}" >> "${PROJDATADIR}/${PROJECT}_list.txt"
					fi
				fi
			fi
		else
			echo "${file} is not a FASTQ(.gz) read file, not acting on it"
		fi
	done
fi

# Task: Invert list so that the important isolates (for us at least) get run first
write_Progress
run_task_id=2
if [[ -f "${PROJDATADIR}/${PROJECT}_list.txt" ]]; then
	sort -k2,2 -t'/' -r "${PROJDATADIR}/${PROJECT}_list.txt" -o "${PROJDATADIR}/${PROJECT}_list.txt"
fi

# Task: Loops through list file to create an array of all isolates to run through pipeline
write_Progress
run_task_id=3
declare -a isolate_list=()
while IFS= read -r file || [ -n "$file" ]; do
	echo "Found: ${file}"
	file=$(echo "${file}" | tr -d '[:space:]')
	isolate_list+=("${file}")
done < "${list_path}"

# Task: Displays number and names of files found to analyze
write_Progress
run_task_id=4
if [[ ${#isolate_list[@]} -gt 1 ]]; then
	echo "Will analyze these ${#isolate_list[@]} files: ${isolate_list[*]}"
elif [[ ${#isolate_list[@]} -eq 1 ]]; then
	echo "Will analyze this file: ${isolate_list[0]}"
else
	echo "No files found in ${list_path}"
fi
isolate_count=${#isolate_list[@]}
run_start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

#Each file in the list is checked individually for successful completion and added then added to the log for the run
log_dir="${PROJDATADIR}"
log_file="${PROJDATADIR}/${PROJECT}_on_${run_start_time}.log"
command_log_file="${PROJDATADIR}/${PROJECT}_on_${run_start_time}_command.log"
echo -e "Below tools and version number of all singularity calls are printed with the exact command used to run analysis:\n" > "${command_log_file}"
exec > >(tee -a ${log_file}) 2>&1

# Task: Get the time the run started to use as the identifier
write_Progress
run_task_id=5
outarray=()
echo "Run started at ${run_start_time}; Log saved to ${log_file}"
echo "Run started at ${run_start_time}" #> "${log_file}"
outarray+=("${PROJECT} started at ${run_start_time} and saved to ${log_file}")

run_task_id=6
loop_inc=0
for isolate in "${isolate_list[@]}"; do
	write_Progress
	isolate_number=${loop_inc}
	#Time tracker to gauge time used by each step
	totaltime=0
	start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

	# Set arguments to isolate_name PROJECT (miseq run id) and SAMPDATADIR(${output_dir}/PROJECT/isolate_name)
	isolate_name=$(echo "${isolate}" | awk -F/ '{print $2}' | tr -d '[:space:]')
	SAMPDATADIR="${PROJDATADIR}/${isolate_name}"

	echo -e "${isolate} started at ${start_time}\n" #>> "${log_file}"
	echo -e "${version_type}:${version_num}" > "${command_log_file}"
	echo -e "-:-:-:-:-:   Starting on isolate: ${isolate_name} :-:-:-:-:-\n" >> "${command_log_file}"


	# Remove old run stats as the presence of the file indicates run completion
	if [[ -f "${SAMPDATADIR}/${isolate_name}_pipeline_stats.txt" ]]; then
		rm "${SAMPDATADIR}/${isolate_name}_pipeline_stats.txt"
	fi

	# Create an empty time_summary file that tracks clock time of tools used
	touch "${SAMPDATADIR}/${isolate_name}_time_summary.txt"
	time_summary="${SAMPDATADIR}/${isolate_name}_time_summary.txt"

	echo "Time summary for ${PROJECT}/${isolate_name}: Started ${global_time}" >> "${time_summary}"
	echo "${PROJECT}/${isolate_name} started at ${global_time}"

	echo "Starting processing of ${PROJECT}/${isolate_name}"
	if [[ "${assemblies}" == "false" ]]; then
		# Task: Checks and prepares FASTQ folder exists for current sample
		write_Progress
		task_number=1
		if [[ -d "${SAMPDATADIR}/FASTQs" ]]; then
			# Checks if FASTQ folder contains any files then continue
			if [[ "$(ls -A "${SAMPDATADIR}/FASTQs")" ]]; then
				# Checks to see if those files in the folder are unzipped fastqs
				count_unzip=`ls -1 ${SAMPDATADIR}/FASTQs/*.fastq 2>/dev/null | wc -l`
				count_zip=`ls -1 ${SAMPDATADIR}/FASTQs/*.fastq.gz 2>/dev/null | wc -l`
				if [[ ${count_unzip} != 0 ]]; then
					echo "----- FASTQ(s) exist, continuing analysis -----"
					if [[ -f "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq" ]] && [[ ! -f "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq.gz" ]]; then
						gzip < "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq" > "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq.gz"
					fi
					if [[ -f "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq" ]] && [[ ! -f "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq.gz" ]]; then
						gzip < "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq" > "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq.gz"
					fi
				# Checks if they are zipped fastqs (checks for R1 first)
				elif [[ -f "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq.gz" ]]; then
					#echo "R1 zipped exists - unzipping"
					gunzip -c "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq.gz" > "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq"
					# Checks for paired R2 file
					if [[ -f "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq.gz" ]]; then
						#echo "R2 zipped exists - unzipping"
						gunzip -c "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq.gz" > "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq"
					else
						echo "No matching R2 to unzip :("
					fi
				# Checks to see if there is an abandoned R2 zipped fastq
				elif [[ -f "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq.gz" ]]; then
					#echo "R2 zipped  exists - unzipping"
					gunzip -c "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq.gz" > "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq"
					echo "No matching R1 to unzip :("
				fi
			# If the folder is empty then return from function
			else
				echo "FASTQs folder empty - No fastqs available for ${isolate_name} (and download was not requested). Either unzip fastqs to ${SAMPDATADIR}/FASTQs or run the -d flag to trigger unzipping of gzs"
			fi
		# If the fastq folder does not exist then return out of function
		else
			echo "FASTQs not downloaded and FASTQs folder does not exist for ${isolate_name}. No fastqs available (and download was not requested). Unzip fastqs to ${SAMPDATADIR}/FASTQs"
		fi

		# Get start time for qc check
		start=$SECONDS
		### Count the number of Q20, Q30, bases and reads within a pair of FASTQ files
		echo "----- Counting read quality -----"
		# Task: Checks for and creates the specified output folder for the QC counts
		write_Progress
		task_number=2
		if [ ! -d "${SAMPDATADIR}/preQCcounts" ]; then
			echo "Creating ${SAMPDATADIR}/preQCcounts"
			mkdir -p "${SAMPDATADIR}/preQCcounts"
		fi
		# Run qc count check on raw reads
		echo -e "${SAMPDATADIR}/${isolate_name}	Q20_Total_[bp]	Q30_Total_[bp]	Q20_R1_[bp]	Q20_R2_[bp]	Q20_R1_[%]	Q20_R2_[%]	Q30_R1_[bp]	Q30_R2_[bp]	Q30_R1_[%]	Q30_R2_[%]	Total_Sequenced_[bp]	Total_Sequenced_[reads]" > "${SAMPDATADIR}/preQCcounts/${isolate_name}_counts.txt"
		python3 "${src}/Fastq_Quality_Printer.py" -1 "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq" -2 "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq" >> "${SAMPDATADIR}/preQCcounts/${isolate_name}_counts.txt"
		# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
		#raw_length_R1=$(cat ${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R1.fsq | paste - - - - | cut -f2 |tr -d '\n' | wc -c)
		#raw_length_R2=$(cat ${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R2.fsq | paste - - - - | cut -f2 |tr -d '\n' | wc -c)
		end=$SECONDS
		timeQCcount=$((end - start))
		echo "QC count - ${timeQCcount} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeQCcount))

		# Task: Trimming and Quality Control
		write_Progress
		task_number=3
		echo "----- Running BBDUK on reads -----"
		# Gets start time for bbduk
		start=$SECONDS
		# Creates folder for BBDUK output
		if [ ! -d "${SAMPDATADIR}/removedAdapters" ]; then
			echo "Creating ${SAMPDATADIR}/removedAdapters"
			mkdir -p "${SAMPDATADIR}/removedAdapters"
		# It complains if a folder already exists, so the current one is removed (shouldnt happen anymore as each analysis will move old runs to new folder)
		else
			echo "Removing old ${SAMPDATADIR}/removedAdapters"
			rm -r "${SAMPDATADIR}/removedAdapters"
			echo "Recreating ${SAMPDATADIR}/removedAdapters"
			mkdir -p "${SAMPDATADIR}/removedAdapters"
		fi
		### Run bbduk
		#singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES ${local_DBs}/singularities/bbtools.simg bbduk.sh -${bbduk_mem} threads=${procs} in=/SAMPDIR/FASTQs/${isolate_name}_R1_001.fastq in2=/SAMPDIR/FASTQs/${isolate_name}_R2_001.fastq out=/SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R1.fsq out2=/SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R2.fsq ref=/DATABASES/phiX.fasta k=${bbduk_k} hdist=${bbduk_hdist}
		#echo -e "bbtools:37.87 -- bbduk.sh -${bbduk_mem} threads=${procs} in=${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq in2=${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq out=${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R1.fsq out2=${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R2.fsq ref=/DATABASES/phiX.fasta k=${bbduk_k} hdist=${bbduk_hdist}" >> "${command_log_file}"
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/staphb/bbtools:38.94 bbduk.sh -${bbduk_mem} threads=${procs} in=/SAMPDIR/FASTQs/${isolate_name}_R1_001.fastq in2=/SAMPDIR/FASTQs/${isolate_name}_R2_001.fastq out=/SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R1.fsq out2=/SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R2.fsq ref=/DATABASES/phiX.fasta k=${bbduk_k} hdist=${bbduk_hdist}
		echo -e "bbtools:39.84 -- bbduk.sh -${bbduk_mem} threads=${procs} in=${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq in2=${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq out=${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R1.fsq out2=${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R2.fsq ref=/DATABASES/phiX.fasta k=${bbduk_k} hdist=${bbduk_hdist}" >> "${command_log_file}"
		# Get end time of bbduk and calculate run time and append to time summary (and sum to total time used)
		remAdapt_length_R1=$(cat ${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R1.fsq | paste - - - - | cut -f2 |tr -d '\n' | wc -c)
		remAdapt_length_R2=$(cat ${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R2.fsq | paste - - - - | cut -f2 |tr -d '\n' | wc -c)
		echo -e "R1:	${remAdapt_length_R1}\nR2:	${remAdapt_length_R2}" > ${SAMPDATADIR}/removedAdapters/no_PhiX_total_lengths.txt


		end=$SECONDS
		timeAdapt=$((end - start))
		echo "Removing Adapters - ${timeAdapt} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeAdapt))


		# Task: Quality and Adapter Trimming using trimmomatic
		write_Progress
		task_number=4
		echo "----- Running fastp on reads -----"
		# Get start time of trimmomatic
		start=$SECONDS
		# Creates folder for trimmomatic output if it does not exist
		if [ ! -d "${SAMPDATADIR}/trimmed" ]; then
			mkdir -p "${SAMPDATADIR}/trimmed"
		fi
		### Run trimmomatic
		#singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/trimmomatic:0.39--1 trimmomatic PE -${trim_phred} -threads ${procs} /SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R1.fsq /SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R2.fsq /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq /SAMPDIR/trimmed/${isolate_name}_R1_001.unpaired.fq /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq /SAMPDIR/trimmed/${isolate_name}_R2_001.unpaired.fq #ILLUMINACLIP:/DATABASES/adapters.fasta:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome} SLIDINGWINDOW:${trim_window_size}:${trim_window_qual} LEADING:${trim_leading} TRAILING:${trim_trailing} MINLEN:${trim_min_length}
		#echo -e "trimmomatic:0.39 -- trimmomatic PE -${trim_phred} -threads ${procs} ${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R1.fsq ${SAMPDATADIR}/removedAdapters/${isolate_name}-noPhiX-R2.fsq ${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq ${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.unpaired.fq ${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq ${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.unpaired.fq ILLUMINACLIP:${local_DBs}/#adapters.fasta:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome} SLIDINGWINDOW:${trim_window_size}:${trim_window_qual} LEADING:${trim_leading} TRAILING:${trim_trailing} MINLEN:${trim_min_length}\n" >> "${command_log_file}"

		### Run fastp instead of trimmomatic
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/fastp:0.23.1--h79da9fb_0 fastp -w ${procs} -i /SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R1.fsq -I /SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R2.fsq -o /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq --unpaired1 /SAMPDIR/trimmed/${isolate_name}.single1.fq -O /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq --unpaired2 /SAMPDIR/trimmed/${isolate_name}.single2.fq --adapter_fasta /DATABASES/adapters.fasta -r --cut_right_window_size ${trim_window_size} --cut_right_mean_quality ${trim_window_qual} -l ${trim_min_length} -g -5 ${trim_leading} -3 ${trim_trailing}
		echo -e "fastp:0.23.1 -- fastp -w ${procs} -i /SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R1.fsq -I /SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R2.fsq -o /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq --unpaired1 /SAMPDIR/trimmed/${isolate_name}.single1.fq -O /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq --unpaired2 /SAMPDIR/trimmed/${isolate_name}.single2.fq --adapter_fasta /DATABASES/adapters.fasta -r --cut_right_window_size ${trim_window_size} --cut_right_mean_quality ${trim_window_qual} -l ${trim_min_length} -g -5 ${trim_leading} -3 ${trim_trailing}"

		# Merge both unpaired fq files into one for GOTTCHA
		cat ${SAMPDATADIR}/trimmed/${isolate_name}.single*.fq > ${SAMPDATADIR}/trimmed/${isolate_name}.single.fq
		rm ${SAMPDATADIR}/trimmed/${isolate_name}.single1.fq
		rm ${SAMPDATADIR}/trimmed/${isolate_name}.single2.fq

		# Get end time of trimmomatic and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeTrim=$((end - start))
		echo "Trimming - ${timeTrim} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeTrim))
		gzip -c "${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq" > "${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq.gz"
		gzip -c "${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq" > "${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq.gz"

		# Task: Check differences after QC and trimming (also for gottcha proper read count for assessing unclassified reads)
		write_Progress
		task_number=5
		# Get start time for qc check on trimmed reads
		start=$SECONDS
		### Count the number of Q20, Q30, bases and reads within the trimmed pair of FASTQ files
		echo "----- Counting read quality of trimmed files-----"
		# Run qc count check on filtered reads
		echo -e "${SAMPDATADIR}/${isolate_name}	Q20_Total_[bp]	Q30_Total_[bp]	Q20_R1_[bp]	Q20_R2_[bp]	Q20_R1_[%]	Q20_R2_[%]	Q30_R1_[bp]	Q30_R2_[bp]	Q30_R1_[%]	Q30_R2_[%]	Total_Sequenced_[bp]	Total_Sequenced_[reads]" > "${SAMPDATADIR}/preQCcounts/${isolate_name}_trimmed_counts.txt"
		python3 "${src}/Fastq_Quality_Printer.py" -1 "${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq" -2 "${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq" >> "${SAMPDATADIR}/preQCcounts/${isolate_name}_trimmed_counts.txt"

		cat "${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq" "${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq" > "${SAMPDATADIR}/trimmed/${isolate_name}.paired.fq"

		# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeQCcount_trimmed=$((end - start))
		echo "QC count trimmed - ${timeQCcount_trimmed} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeQCcount))

		# Task: Run Kraken on cleaned reads  ######
		write_Progress
		task_number=6
		echo "----- Running Kraken on cleaned reads -----"
		# Get start time of kraken on reads
		start=$SECONDS
		# Create directory for kraken dataset
		if [[ ! -d "${SAMPDATADIR}/kraken/preAssembly" ]]; then
			mkdir -p "${SAMPDATADIR}/kraken/preAssembly"
		fi
		# Run kraken
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken --paired --db /DATABASES/kraken/${kraken_DB}  --preload --fastq-input --threads 4 --output /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.kraken --classified-out /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.classified /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq
		echo -e "kraken:1.1.1 -- kraken --paired --db ${local_DBs}/${kraken_DB}  --preload --fastq-input --threads 4 --output ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.kraken --classified-out ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.classified ${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq ${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq\n" >> "${command_log_file}"
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken-mpa-report --db /DATABASES/kraken/${kraken_DB} /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.kraken > "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.mpa"
		echo -e "kraken:1.1.1 -- kraken-mpa-report --db ${local_DBs}/${kraken_DB} ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.kraken > ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.mpa\n" >> "${command_log_file}"
		python3 "${src}/Metaphlan2krona.py" -p "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.mpa" -k "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.krona"
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR docker://quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2 ktImportText /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.krona -o /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.html
		echo -e"krona:2.8 -- ktImportText ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.krona -o ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.html\n" >> "${command_log_file}"
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken-report --db /DATABASES/kraken/${kraken_DB} /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.kraken > "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.list"
		echo -e "kraken:1.1.1 -- kraken-report --db ${local_DBs}/${kraken_DB} ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.kraken > ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.list\n" >> "${command_log_file}"
		"${src}/best_hit_from_kraken.sh" "${SAMPDATADIR}" "pre" "paired" "kraken"
		# Get end time of kraken on reads and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeKrak=$((end - start))
		echo "Kraken - ${timeKrak} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeKrak))

		# # Task: Run gottcha(v1) on cleaned reads
		# write_Progress
		# task_number=7
		# echo "----- Running gottcha on cleaned reads -----"
		# # Get start time of gottcha
		# start=$SECONDS
		# # run gottcha
		# singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B "${local_DBs}":/DATABASES ${local_DBs}/singularities/gottcha.simg gottcha.pl --mode all --outdir /SAMPDIR/gottcha/gottcha_S --input /SAMPDIR/trimmed/${isolate_name}.paired.fq --database /DATABASES/gottcha/${gottcha_DB}
		# echo -e "gottcha:1.0b -- gottcha.pl --mode all --outdir ${SAMPDATADIR}/gottcha/gottcha_S --input ${SAMPDATADIR}/trimmed/${isolate_name}.paired.fq --database ${local_DBs}/gottcha/${gottcha_DB}\n" >> "${command_log_file}"
		# ### Public version failed BWA
		# #singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B "${local_DBs}":/DATABASES docker://quay.io/biocontainers/gottcha:1.0--pl526_2 gottcha.pl --mode all --outdir /SAMPDIR/gottcha/gottcha_S --input /SAMPDIR/trimmed/${isolate_name}.paired.fq --database /DATABASES/gottcha/${gottcha_DB}
		# #echo -e "gottcha:1.0b -- gottcha.pl --mode all --outdir ${SAMPDATADIR}/gottcha/gottcha_S --input ${SAMPDATADIR}/trimmed/${isolate_name}.paired.fq --database ${local_DBs}/gottcha/${gottcha_DB}\n" >> "${command_log_file}"
		#
		# singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR docker://quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2 ktImportText /SAMPDIR/gottcha/gottcha_S/${isolate_name}_temp/${isolate_name}.lineage.tsv -o /SAMPDIR/gottcha/${isolate_name}_species.krona.html
		# echo -e "krona:2.8 -- ktImportText ${SAMPDATADIR}/gottcha/gottcha_S/${isolate_name}_temp/${isolate_name}.lineage.tsv -o ${SAMPDATADIR}/gottcha/${isolate_name}_species.krona.html\n" >> "${command_log_file}"
		#
		# #Create a best hit from gottcha1 file
		# "${src}/best_hit_from_gottcha1.sh" "${SAMPDATADIR}"
		# # Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
		# end=$SECONDS
		# timeGott=$((end - start))
		# echo "Gottcha - ${timeGott} seconds" >> "${time_summary}"
		# totaltime=$((totaltime + timeGott))

		# Task: Check reads using SRST2
		write_Progress
		task_number=7
		echo "----- Running SRST2 -----"
		start=$SECONDS

		if [[ ! -d "${SAMPDATADIR}/srst2" ]]; then
				mkdir "${SAMPDATADIR}/srst2"
		fi

		# prep file for srst transfer
		cp ${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz
		cp ${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz

    # create temp location for srst2 db
    mkdir -p "${SAMPDATADIR}/srst2/tmp_db"
    cp "${local_DBs}/star/${ResGANNCBI_srst2_filename}_srst2.fasta" "${SAMPDATADIR}/srst2/tmp_db/"

		singularity exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DBs ${local_DBs}/singularities/srst2.simg python /usr/local/bin/srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/${isolate_name} --gene_db /SAMPDIR/srst2/tmp_db/${ResGANNCBI_srst2_filename}_srst2.fasta
		echo -e "srst2:0.2.0 -- srst2 --input_pe ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output ${SAMPDATADIR}/srst2/${isolate_name}_ResGANNCBI --gene_db ${SAMPDATADIR}/srst2/tmp_db/${ResGANNCBI_srst2_filename}_srst2.fasta\n" >> "${command_log_file}"

		# Cleans up leftover files
		rm "${SAMPDATADIR}/srst2/"*".bam"
		rm "${SAMPDATADIR}/srst2/"*".pileup"
		rm "${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz"
		rm "${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz"
    rm -r "${SAMPDATADIR}/srst2/tmp_db"

		# # Removes the extra ResGANNCBI__ from all files created
		# find ${SAMPDATADIR}/srst2 -type f -name "*ResGANNCBI__*" | while read FILE ; do
		#   dirname=$(dirname $FILE)
		# 	filename=$(basename $FILE)
		# 	filename="${filename/_ResGANNCBI__/__}"
		# 	#echo "Found-${FILE}"
		# 	#echo "${filename}"
		#     mv "${FILE}" "${dirname}/${filename}"
		# done

		end=$SECONDS
		timesrst2=$((end - start))
		echo "SRST2 - ${timesrst2} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timesrst2))

		# Task: Assembling Using SPAdes
		write_Progress
		task_number=8
		echo "----- Assembling Using SPAdes -----"
		# Get start time of SPAdes
		start=$SECONDS
		# script tries 3 times for a completed assembly
		for i in 1 2 3
		do
			# If assembly exists already and this is the first attempt (then the previous run will be used) [should not happen anymore as old runs are now renamed]
			if [ -s "${SAMPDATADIR}/Assembly/scaffolds.fasta" ]; then
				echo "Previous assembly already exists, using it (delete/rename the assembly folder at ${SAMPDATADIR}/ if you'd like to try to reassemble"
			# Run normal mode if no assembly file was found
			else
				singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/spades:3.15.3--h95f258a_0 spades.py --careful --only-assembler --pe1-1 /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq --pe1-2 /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq --pe1-s /SAMPDIR/trimmed/${isolate_name}.single.fq -o /SAMPDIR/Assembly --phred-offset "${spades_phred_offset}" -t "${procs}"
				echo -e "SPAdes:3.15.3 -- spades.py --careful --only-assembler --pe1-1 ${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq --pe1-2 ${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq --pe1-s ${SAMPDATADIR}/trimmed/${isolate_name}.single.fq -o ${SAMPDATADIR}/Assembly --phred-offset ${spades_phred_offset} -t ${procs}\n" >> "${command_log_file}"
			fi
			# Removes any core dump files (Occured often during testing and tweaking of memory parameter
			if [ -n "$(find "${src}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
				echo "Found core dump files in assembly (assumed to be from SPAdes, but could be anything before that as well) and attempting to delete"
				find "${src}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
			fi
		done
		# Returns if all 3 assembly attempts fail
		if [[ -f "${SAMPDATADIR}/Assembly/scaffolds.fasta" ]] && [[ -s "${SAMPDATADIR}/Assembly/scaffolds.fasta" ]]; then
			echo "Assembly completed and created a non-empty scaffolds file"
		else
			echo "Assembly FAILED 3 times, continuing to next sample..." >&2
			return 1
		fi
		# Get end time of SPAdes and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeSPAdes=$((end - start))
		echo "SPAdes - ${timeSPAdes} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeSPAdes))

		# Task: Removing Short Contigs
		write_Progress
		task_number=9
		echo "----- Removing Short Contigs -----"
		python3 "${src}/removeShortContigs.py" -i "${SAMPDATADIR}/Assembly/scaffolds.fasta" -t 500 -s "normal_SPAdes"
		mv "${SAMPDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta"
	fi

	# Resumes remaineder of pipeline if assemblies were set to true or retry
	# Checks to see that the trimming and renaming worked properly, returns if unsuccessful
	if [ ! -s "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" ]; then
		echo "Trimmed contigs file does not exist continuing to next sample">&2
		continue
	fi

	# Task: ReKraken on Assembly
	write_Progress
	task_number=10
	echo "----- Running Kraken on Assembly -----"
	# Get start time of kraken on assembly
	start=$SECONDS
	# Create directory for kraken dataset on assembly
	if [[ ! -d "${SAMPDATADIR}/kraken/postAssembly" ]]; then
		mkdir -p "${SAMPDATADIR}/kraken/postAssembly"
	fi
	# Run kraken on assembly
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken --db /DATABASES/kraken/${kraken_DB}  --preload --threads ${procs} --output /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.kraken --classified-out /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.classified /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta
	echo -e "kraken:1.1.1 -- kraken --db ${local_DBs}/${kraken_DB}  --preload --threads ${procs} --output ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.kraken --classified-out ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.classified ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta\n"  >> "${command_log_file}"
	python3 ${src}/Kraken_Assembly_Converter_2_Exe.py -i "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.kraken"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken-mpa-report --db /DATABASES/kraken/${kraken_DB} /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled_BP.kraken > ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted.mpa
  echo -e "kraken:1.1.1 -- kraken-mpa-report --db ${local_DBs}/${kraken_DB} ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_BP.kraken > ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted.mpa\n" >> "${command_log_file}"
	python3 "${src}/Metaphlan2krona.py" -p "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted.mpa" -k "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted.krona"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken-report --db /DATABASES/kraken/${kraken_DB} /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled_BP.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_BP.list"
	echo -e "kraken:1.1.1 -- kraken-report --db ${local_DBs}/${kraken_DB} ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_BP.kraken > ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_BP.list\n" >> "${command_log_file}"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2 ktImportText /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled_weighted.krona -o /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled_weighted_BP_krona.html
	echo -e "krona:2.8 -- ktImportText ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted.krona -o ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted_BP_krona.html\n" >> "${command_log_file}"
	"${src}/best_hit_from_kraken.sh" "${SAMPDATADIR}" "post" "assembled_BP" "kraken"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken-report --db /DATABASES/kraken/${kraken_DB} /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.list"
	echo -e "kraken:1.1.1 -- kraken-report --db ${local_DBs}/${kraken_DB} ${SAMPDATADIR}/${isolate_name}_assembled.kraken > ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.list\n" >> "${command_log_file}"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken-mpa-report --db /DATABASES/kraken/${kraken_DB} /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.mpa"
	echo -e "kraken:1.1.1 -- kraken-mpa-report --db ${local_DBs}/${kraken_DB} ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.kraken > ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.mpa\n" >> "${command_log_file}"
	python3 "${src}/Metaphlan2krona.py" -p "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.mpa" -k "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.krona"
	singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR docker://quay.io/biocontainers/krona:2.8--pl5262hdfd78af_2 ktImportText /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.krona -o /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.html
	echo -e "krona:2.8 -- ktImportText ${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.krona -o ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.html\n" >> "${command_log_file}"
	singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.1.1--pl5262h7d875b9_5 kraken-report --db /DATABASES/kraken/${kraken_DB} /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.list"
	echo -e "kraken:1.1.1 -- kraken-report --db ${local_DBs}/${kraken_DB} ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.kraken > ${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.list\n" >> "${command_log_file}"
	"${src}/best_hit_from_kraken.sh" "${SAMPDATADIR}" "post" "assembled" "kraken"
	# Get end time of kraken on assembly and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeKrakAss=$((end - start))
	echo "Kraken on Assembly - ${timeKrakAss} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeKrakAss))

	# Task:  Get ID fom 16s
	write_Progress
	task_number=11
	echo "----- Identifying via 16s blast -----"
	start=$SECONDS

	if [ ! -d "${SAMPDATADIR}/16s" ]; then
		echo "Creating ${SAMPDATADIR}/16s"
		mkdir "${SAMPDATADIR}/16s"
	fi

	# Get original working directory so that it can return to it after running (may not need to do this but havent tested it out yet)
	owd=$(pwd)
	cd ${SAMPDATADIR}/16s

	# Run barrnap to discover ribosomal sequences
	#barrnap --kingdom bac --threads ${procs} "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" > ${SAMPDATADIR}/16s/${isolate_name}_rRNA_finds.txt

	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/barrnap:0.9--hdfd78af_4 barrnap --kingdom bac /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta > "${SAMPDATADIR}/16s/${isolate_name}_rRNA_finds.txt"
	echo -e "barrnap:0.9 -- barrnap --kingdom bac ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta > ${SAMPDATADIR}/16s/${isolate_name}_rRNA_finds.txt\n"  >> "${command_log_file}"
	# Checks for successful output from barrnap, *rRNA_seqs.fasta
	if [[ ! -s ${SAMPDATADIR}/16s/${isolate_name}_rRNA_finds.txt ]]; then
		echo "rNA_seqs.fasta does NOT exist"
	fi

	# Checks barrnap output and finds all 16s hits and creates a multi-fasta file to list all possible matches
	lines=0
	found_16s="false"
	while IFS='' read -r line || [ -n "$line" ]; do
		contig=$(echo ${line} | cut -d' ' -f1)
		cstart=$(echo ${line} | cut -d' ' -f4)
		cstop=$(echo ${line} | cut -d' ' -f5)
		ribosome=$(echo ${line} | cut -d' ' -f9 | cut -d'=' -f3)
		if [ "${ribosome}" = "16S" ]; then
			# Replace with subsequence once it can handle multi-fastas
			#make_fasta $1 $2 $contig $cstart $cstop
			python3 ${src}/get_subsequence.py -i "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" -s ${cstart} -e ${cstop} -t ${contig} -o "${contig} 16s-${lines}" >> ${SAMPDATADIR}/16s/${isolate_name}_16s_rna_seq_${lines}.fasta
			found_16s="true"
			lines=$((lines + 1))
		fi
	done < "${SAMPDATADIR}/16s/${isolate_name}_rRNA_finds.txt"

	# Adds No hits found to output file in the case where no 16s ribosomal sequences were found
	if [[ "${found_16s}" == "false" ]]; then
		echo -e "best_hit	${isolate_name}	No_16s_sequences_found" > "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt"
		echo -e "largest_hit	${isolate_name}	No_16s_sequences_found" >> "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt"
	fi

	# Blasts the NCBI database to find the closest hit to every entry in the 16s fasta list
	###### MAX_TARGET_SEQS POSSIBLE ERROR
	line_inc=0
	while [[ -f ${SAMPDATADIR}/16s/${isolate_name}_16s_rna_seq_${line_inc}.fasta ]]; do
		echo "Blasting ${isolate_name}_16s_rna_seq_${line_inc}.fasta"
		singularity exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0 blastn -word_size 10 -task blastn -remote -db nt -max_hsps 1 -max_target_seqs 10 -query /SAMPDIR/16s/${isolate_name}_16s_rna_seq_${line_inc}.fasta -out /SAMPDIR/16s/${isolate_name}.nt.RemoteBLASTN_${line_inc} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ssciname"
		echo "blast:2.12.0 blastn -word_size 10 -task blastn -remote -db nt -max_hsps 1 -max_target_seqs 1 -query /SAMPDIR/16s/${isolate_name}_16s_rna_seq_${line_inc}.fasta -out /SAMPDIR/16s/${isolate_name}.nt.RemoteBLASTN_${line_inc} -outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ssciname"
		line_inc=$(( line_inc + 1 ))
	done

	cat ${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_* > ${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_all

	# Sorts the list based on sequence match length to find the largest hit
	sort -k4 -n "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_all" --reverse > "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_all.sorted"

	# Gets taxon info from the best bitscore (literal top) hit from the blast list
	if [[ -s "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_all" ]]; then
		me=$(whoami)
		accessions=$(head -n 1 "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_all")
		hits=$(echo "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_all" | wc -l)
	#	echo ${accessions}
		gb_acc=$(echo "${accessions}" | cut -d'|' -f4)
		echo ${gb_acc}
		attempts=0
		# Will try getting info from entrez up to 5 times, as it has a higher chance of not finishing correctly on the first try
		while [[ ${attempts} -lt 5 ]]; do
			blast_id=$(python ${src}/entrez_get_taxon_from_accession.py -a "${gb_acc}" -e "${me}@cdc.gov")
			if [[ ! -z ${blast_id} ]]; then
				break
			else
				attempts=$(( attempts + 1 ))
			fi
			sleep 1
		done
		echo ${blast_id}
		if [[ -z ${blast_id} ]]; then
			blast_id="No_16s_matches_found"
		fi
		#blast_id=$(echo ${blast_id} | tr -d '\n')
		echo -e "best_hit	${isolate_name}	${blast_id}" > "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt"
		if [[ "${hits}" -eq 1 ]]; then
			echo -e "largest	${isolate_name}	${blast_id}" >> "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt"
			skip_largest="true"
		fi
	else
		echo "No remote blast file"
	fi

	best_blast_id=${blast_id}

	# Gets taxon info from the largest hit from the blast list
	if [[ ${skip_largest} != "true" ]]; then
		# Gets taxon info from the largest hit from the blast list
		if [[ -s "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_all.sorted" ]]; then
			me=$(whoami)
			accessions=$(head -n 1 "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN_all.sorted")
			gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
			attempts=0
			# Will try getting info from entrez up to 5 times, as it has a higher chance of not finishing correctly on the first try
			while [[ ${attempts} -lt 5 ]]; do
				blast_id=$(python ${shareScript}/entrez_get_taxon_from_accession.py -a "${gb_acc}" -e "${me}@cdc.gov")
				if [[ ! -z ${blast_id} ]]; then
					break
				else
					attempts=$(( attempts + 1 ))
				fi
			done
			echo ${blast_id}
			if [[ -z ${blast_id} ]]; then
				blast_id="No_16s_matches_found"
			fi
			#	blast_id$(echo ${blast_id} | tr -d '\n')
			if [[ "${hits}" -eq 1 ]] && [[ "${best_blast_id}" == "No_16s_matches_found" ]]; then
				echo -e "best_hit	${isolate_name}	${blast_id}" > "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt"
			fi
			echo -e "largest	${isolate_name}	${blast_id}" >> "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt"
		else
			echo "No sorted remote blast file"
		fi
	fi

	# Go back to original working directory
	cd ${owd}
	end=$SECONDS
	time16s=$((end - start))
	echo "16S - ${time16s} seconds" >> "${time_summary}"
	totaltime=$((totaltime + time16s))

	# Task: Check quality of Assembly
	write_Progress
	task_number=12
	echo "----- Running quality checks on Assembly -----"
	# Get start time of QC assembly check
	start=$SECONDS
	# Run qc assembly check

	owd="$(pwd)"
 # Checks for output folder existence and creates creates if not
	if [ ! -d "${SAMPDATADIR}/Assembly_Stats" ]; then
		echo "Creating ${SAMPDATADIR}/Assembly_Stats"
		mkdir -p "${SAMPDATADIR}/Assembly_Stats"
	fi
	cd "${SAMPDATADIR}/Assembly_Stats"
	# Call QUAST
	#python2 "/apps/x86_64/quast/quast-4.3/quast.py" -o "${SAMPDATADIR}/Assembly_Stats" "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta"
	#singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${local_DBs}/singularities/QUAST5.simg python3 /quast/quast.py --no-icarus --no-html --no-snps -o /SAMPDIR/Assembly_Stats /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta
	#echo -e "Quast:5.0.0 -- python3 /quast/quast.py -o ${SAMPDATADIR}/Assembly_Stats ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta\n" >> "${command_log_file}"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/quast:5.0.2--py37pl5262hfecc14a_5 quast.py --no-icarus --no-html --no-snps -o /SAMPDIR/Assembly_Stats /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta
	echo -e "Quast:5.0.2 -- python3 quast.py -o ${SAMPDATADIR}/Assembly_Stats ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta\n" >> "${command_log_file}"
	mv "${SAMPDATADIR}/Assembly_Stats/report.txt" "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.txt"
	mv "${SAMPDATADIR}/Assembly_Stats/report.tsv" "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv"
	# Get end time of qc quality check and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeQCcheck=$((end - start))
	echo "QC Check of Assembly - ${timeQCcheck} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeQCcheck))
	cd "${owd}"

	# Task: Prokka on assembly
	write_Progress
	task_number=13
	echo "----- Running Prokka on Assembly -----"
	# Get start time for prokka
	start=$SECONDS
	# Run prokka
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/prokka:1.14.6--pl5262hdfd78af_1 prokka --outdir /SAMPDIR/prokka /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta
	echo -e "prokka:1.14.6 -- prokka --outdir ${SAMPDATADIR}/prokka ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta\n" >> "${command_log_file}"
	#echo "About to rename files"
	for pfile in ${SAMPDATADIR}/prokka/*.*; do
		fullname=$(basename "${pfile}")
		ext="${fullname##*.}"
		echo "Renaming ${pfile} to ${SAMPDATADIR}/prokka/${isolate_name}_PROKKA.${ext}"
		mv "${pfile}" "${SAMPDATADIR}/prokka/${isolate_name}_PROKKA.${ext}"
	done
	# Get end time of prokka and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeProk=$((end - start))
	echo "Identifying Genes (Prokka) - ${timeProk} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeProk))

	# Task:  Rename contigs to something helpful (Had to wait until after prokka runs due to the strict naming requirements
	write_Progress
	task_number=14
	mv "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed_original.fasta"
	python3 "${src}/fasta_headers.py" -i "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed_original.fasta" -o "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta"

	# Task: Average Nucleotide Identity
	write_Progress
	task_number=15
	echo "----- Running ANI for Species confirmation -----"
	# ANI uses assembly and sample would have exited already if assembly did not complete, so no need to check
	# Get start time of ANI
	start=$SECONDS
	# run ANI

	# Checks to see if the local DB ANI folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
	if [ ! -d "${SAMPDATADIR}/ANI/localANIDB_REFSEQ" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/ANI/localANIDB_REFSEQ"
		mkdir -p "${SAMPDATADIR}/ANI/localANIDB_REFSEQ"
	else
		rm -r "${SAMPDATADIR}/ANI/localANIDB_REFSEQ"
		mkdir -p "${SAMPDATADIR}/ANI/localANIDB_REFSEQ"
	fi

	#Copies the samples assembly contigs to the local ANI db folder
	cp "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" "${SAMPDATADIR}/ANI/localANIDB_REFSEQ/sample.fasta"

	#mash dist "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" "${REFSEQ}" > "${SAMPDATADIR}/ANI/${isolate_name}_${REFSEQ_date}_mash.dists"
	singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/mash:2.3--he348c14_1 mash dist /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta /DATABASES/ANI/REFSEQ_${REFSEQ_date}_Bacteria_complete.msh > "${SAMPDATADIR}/ANI/${isolate_name}_${REFSEQ_date}_mash.dists"
	echo -e "mash:2.3 -- mash dist ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta ${local_DBs}/ANI/REFSEQ_${REFSEQ_date}_Bacteria_complete.msh > ${SAMPDATADIR}/ANI/${isolate_name}_${REFSEQ_date}_mash.dists\n" >> "${command_log_file}"
	sort -k3 -n -o "${SAMPDATADIR}/ANI/${isolate_name}_${REFSEQ_date}_mash_sorted.dists" "${SAMPDATADIR}/ANI/${isolate_name}_${REFSEQ_date}_mash.dists"
	#rm "${SAMPDATADIR}/ANI/${isolate_name}_${REFSEQ_date}_mash.dists"

	cutoff=$(head -n${max_ani_samples} "${SAMPDATADIR}/ANI/${isolate_name}_${REFSEQ_date}_mash_sorted.dists" | tail -n1 | cut -d'	' -f3)

	echo "Cutoff IS: ${cutoff}"

	while IFS= read -r var; do
		echo "${var}"
		dist=$(echo ${var} | cut -d' ' -f3)
		kmers=$(echo ${var} | cut -d' ' -f5 | cut -d'/' -f1)
		echo "dist-${dist}"
		if (( $(echo "$dist <= $cutoff" | bc -l) )) && [ ${kmers} -gt 0 ]; then
			filename=$(echo ${var} | cut -d' ' -f2 | rev | cut -d'/' -f1 | rev | cut -d'_' -f3- | rev | cut -d'_' -f2,3,4 | rev)
			alpha=${filename:4:3}
			beta=${filename:7:3}
			charlie=${filename:10:3}
			echo "Copying - ${filename}"
			echo "${filename}" >> ${SAMPDATADIR}/ANI/REFSEQ_${REFSEQ_date}_accessions.txt
			echo "Trying - wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz -P ${SAMPDATADIR}/ANI/localANIDB_REFSEQ"
			wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz -P ${SAMPDATADIR}/ANI/localANIDB_REFSEQ
		else
			break
		fi
	done < ${SAMPDATADIR}/ANI/${isolate_name}_${REFSEQ_date}_mash_sorted.dists

	#mv ${src}/*_genomic.fna.gz "${SAMPDATADIR}/ANI/localANIDB_REFSEQ/"

	#"${src}/append_taxonomy_to_ncbi_assembly_filenames.sh" "${SAMPDATADIR}/ANI/localANIDB_REFSEQ"

	for i in ${SAMPDATADIR}/ANI/localANIDB_REFSEQ/*.gz; do
		old_name=$(basename ${i} | cut -d'.' -f1,2)
		new_name=$(echo ${old_name} | tr -d '[],')
		dir_name=$(dirname ${i})
		gunzip ${i}
		tax_genus=$(head -n1 "${dir_name}/${old_name}.fna" | cut -d' ' -f2 | tr -d '[],')
		tax_species=$(head -n1 "${dir_name}/${old_name}.fna" | cut -d' ' -f3 | tr -d '[],')
		echo "Taxes: ${tax_genus}:${tax_species}"
		mv ${dir_name}/${old_name}.fna ${dir_name}/${tax_genus}_${tax_species}_${new_name}.fasta
	done

	# #Renames all files in the localANIDB_REFSEQ folder by changing extension from fna to fasta (which pyani needs)
	# for file in ${SAMPDATADIR}/ANI/localANIDB_REFSEQ/*.fna;
	# do
	# 	fasta_name=$(basename "${file}" .fna)".fasta"
	# 	mv "${file}" "${SAMPDATADIR}/ANI/localANIDB_REFSEQ/${fasta_name}"
	# done

	# Checks for a previous copy of the aniM folder, removes it if found
	if [ -d "${SAMPDATADIR}/ANI/aniM_REFSEQ" ]; then
		echo "Removing old ANIm results in ${SAMPDATADIR}/ANI/aniM_REFSEQ"
		rm -r "${SAMPDATADIR}/ANI/aniM_REFSEQ"
	fi

	#Calls pyani on local db folder
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/pyani:0.2.11--pyhdfd78af_0 average_nucleotide_identity.py -i /SAMPDIR/ANI/localANIDB_REFSEQ -o /SAMPDIR/ANI/aniM_REFSEQ
	echo -e "pyani:0.2.11 -- average_nucleotide_identity.py -i ${SAMPDATADIR}/ANI/localANIDB_REFSEQ -o ${SAMPDATADIR}/ANI/aniM_REFSEQ\n" >> "${command_log_file}"



	#Extracts the query sample info line for percentage identity from the percent identity file
	while IFS='' read -r line; do
	#	echo "!-${line}"
		if [[ ${line:0:6} = "sample" ]]; then
			sample_identity_line=${line}
	#		echo "found it-"$sample_identity_line
			break
		fi
	done < "${SAMPDATADIR}/ANI/aniM_REFSEQ/ANIm_percentage_identity.tab"

	#Extracts the query sample info line for percentage identity from the percent identity file
	while IFS='' read -r line; do
	#	echo "!-${line}"
		if [[ ${line:0:6} = "sample" ]]; then
			sample_coverage_line=${line}
	#		echo "found it-"$sample_identity_line
			break
		fi
	done < "${SAMPDATADIR}/ANI/aniM_REFSEQ/ANIm_alignment_coverage.tab"

	#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
	if [[ -s "${SAMPDATADIR}/ANI/aniM_REFSEQ/ANIm_percentage_identity.tab" ]]; then
		header_line=$(head -n 1 "${SAMPDATADIR}/ANI/aniM_REFSEQ/ANIm_percentage_identity.tab")
	else
		echo "No ${SAMPDATADIR}/ANI/aniM_REFSEQ/ANIm_percentage_identity.tab file, exiting"
		exit 1
	fi

	#Arrays to read sample names and the %ids for the query sample against those other samples
	IFS="	" read -r -a samples <<< "${header_line}"
	IFS="	" read -r -a percents <<< "${sample_identity_line}"
	IFS="	" read -r -a coverages <<< "${sample_coverage_line}"

	#How many samples were compared
	n=${#samples[@]}

	#Extracts all %id against the query sample (excluding itself) and writes them to file
	for (( i=0; i<n; i++ ));
	do
	#	echo ${i}-${samples[i]}
		if [[ ${samples[i]:0:6} = "sample" ]];
		then
	#		echo "Skipping ${i}"
			continue
		fi
		definition=$(head -1 "${SAMPDATADIR}/ANI/localANIDB_REFSEQ/${samples[i]}.fasta")
		# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
		echo "${percents[i+1]}	${coverages[i+1]}	${samples[i]}	${definition}" >> "${SAMPDATADIR}/ANI/best_hits.txt"
	done

	#Sorts the list in the file based on %id (best to worst)
	sort -nr -t' ' -k1 -o "${SAMPDATADIR}/ANI/best_hits_ordered.txt" "${SAMPDATADIR}/ANI/best_hits.txt"
	#Extracts the first line of the file (best hit)
	best=$(head -n 1 "${SAMPDATADIR}/ANI/best_hits_ordered.txt")
	#Creates an array from the best hit
	IFS='	' read -r -a def_array <<< "${best}"
	#echo -${def_array[@]}+
	#Captures the assembly file name that the best hit came from
	best_file=${def_array[2]}
	#Formats the %id to standard percentage (xx.xx%)
	best_percent=$(awk -v per="${def_array[0]}" 'BEGIN{printf "%.2f", per * 100}')
	best_coverage=$(awk -v per="${def_array[1]}" 'BEGIN{printf "%.2f", per * 100}')
	#echo "${best_file}"

	# Pulling taxonomy from filename which was looked up. Can possibly be out of date. REFSEQ file will ALWAYS be current though
	best_genus=$(echo "${best_file}" | cut -d'_' -f1)
	best_species=$(echo "${best_file}" | cut -d'_' -f2)
	best_organism_guess="${best_genus} ${best_species}"

	#Creates a line at the top of the file to show the best match in an easily readable format that matches the style on the MMB_Seq log
	echo -e "${best_percent}%ID-${best_coverage}%COV-${best_organism_guess}(${best_file}.fna)\\n$(cat "${SAMPDATADIR}/ANI/best_hits_ordered.txt")" > "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_REFSEQ_${REFSEQ_date}).txt"

	#Removes the transient hit files
	if [ -s "${SAMPDATADIR}/ANI/best_hits.txt" ]; then
		rm "${SAMPDATADIR}/ANI/best_hits.txt"
	#	echo "1"
	fi
	if [ -s "${SAMPDATADIR}/ANI/best_hits_ordered.txt" ]; then
		rm "${SAMPDATADIR}/ANI/best_hits_ordered.txt"
	#	echo "2"
	fi

	# Get end time of ANI and calculate run time and append to time summary (and sum to total time used
	end=$SECONDS
	timeANI=$((end - start))
	echo "autoANI - ${timeANI} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeANI))

	# Task: Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
	write_Progress
	task_number=16
	"${src}/determine_taxID.sh" "${SAMPDATADIR}" "${local_DBs}"
	# Capture the anticipated taxonomy of the sample using kraken on assembly output
	echo "----- Extracting Taxonomy from Taxon Summary -----"
	# Checks to see if the kraken on assembly completed successfully
	if [ -s "${SAMPDATADIR}/${isolate_name}.tax" ]; then
		# Read each line of the kraken summary file and pull out each level  taxonomic unit and store for use later in busco and ANI
		while IFS= read -r line; do
			# Grab first letter of line (indicating taxonomic level)
			first=${line::1}
			# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
			if [ "${first}" = "s" ]
			then
				species=$(echo "${line}" | awk -F '	' '{print $2}')
			elif [ "${first}" = "G" ]
			then
				genus=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "F" ]
			then
				family=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "O" ]
			then
				order=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "C" ]
			then
				class=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "P" ]
			then
				phylum=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "K" ]
			then
				kingdom=$(echo "${line}" | awk -F ' ' '{print $2}')
			elif [ "${first}" = "D" ]
			then
				domain=$(echo "${line}" | awk -F ' ' '{print $2}')
			fi
		done < "${SAMPDATADIR}/${isolate_name}.tax"
		# Print out taxonomy for confirmation/fun
		echo "Taxonomy - ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species}"
		# If no kraken summary file was found
	else
		echo "No Taxonomy output available to make best call from, skipped"
	fi

	# Task: BUSCO on prokka output
	write_Progress
	task_number=17
	echo "----- Running BUSCO on Assembly -----"
	# Check to see if prokka finished successfully
	if [ -s "${SAMPDATADIR}/prokka/${isolate_name}_PROKKA.gbf" ] || [ -s "${SAMPDATADIR}/prokka/${isolate_name}_PROKKA.gff" ]; then
		# Get start time of busco
		start=$SECONDS
		# Set default busco database as bacteria in event that we dont have a database match for sample lineage
		buscoDB="bacteria_odb10"
		#buscoDB=$(find ${local_DBs}/BUSCO/ -type d -name "bacteria_odb1"*)
		# Iterate through taxon levels (species to domain) and test if a match occurs to entry in database. If so, compare against it
		busco_found=0
		for tax in $species $genus $family $order $class $phylum $kingdom $domain
		do
			# Check for problematic taxonomy naming
			if [ "${tax}" = "actinobacteria" ]; then
				if [ "${class}" = "actinobacteria" ]; then
					tax="actinobacteria_class"
				elif [ "${phylum}" = "actinobacteria" ]; then
					tax="actinobacteria_phylum"
				fi
			fi
			if [ -d "${local_DBs}/BUSCO/${tax,}_odb10" ]; then
				buscoDB="${tax,}_odb10"
				busco_found=1
				break
			fi
		done
		# Report an unknown sample to the maintenance file to look into
		if [[ "${busco_found}" -eq 0 ]]; then
			global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
			#echo "BUSCO: ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species} - Found as ${PROJECT}/${isolate_name} on ${global_time}" >> "${src}/maintenance_To_Do.txt"
		fi
		# Show which database entry will be used for comparison
		echo "buscoDB:${buscoDB}"
		# Run busco
		# Odd and annoying personal fix for BUSCO strange behaviour
		singularity -s exec -B ${src}/busco_config.ini:/usr/local/config/config.ini -B ${SAMPDATADIR}:/output -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}/BUSCO:/DATABASES docker://quay.io/biocontainers/busco:3.0.2--py35_4 run_BUSCO.py -i /SAMPDIR/prokka/${isolate_name}_PROKKA.faa -o ${isolate_name}_BUSCO -l /DATABASES/${buscoDB} -m prot
		echo -e "BUSCO:3.0.2 -- run_BUSCO.py -i ${SAMPDATADIR}/prokka/${isolate_name}_PROKKA.faa -o ${isolate_name}_BUSCO -l ${local_DBs}/BUSCO/${buscoDB} -m prot\n" >> "${command_log_file}"
		if [ ! -d "${SAMPDATADIR}/BUSCO" ]; then  #create outdir if absent
			echo "Creating ${SAMPDATADIR}/BUSCO"
			mkdir -p "${SAMPDATADIR}/BUSCO"
		fi
		mv ${SAMPDATADIR}/run_${isolate_name}_BUSCO/* ${SAMPDATADIR}/BUSCO/
		### Different naming than normal pipeline
		mv ${SAMPDATADIR}/BUSCO/short_summary_${isolate_name}_BUSCO.txt ${SAMPDATADIR}/BUSCO/short_summary_${isolate_name}.txt
		mv ${SAMPDATADIR}/BUSCO/full_table_${isolate_name}_BUSCO.tsv ${SAMPDATADIR}/BUSCO/full_table_${isolate_name}.tsv
		mv ${SAMPDATADIR}/BUSCO/missing_busco_list_${isolate_name}_BUSCO.tsv ${SAMPDATADIR}/BUSCO/missing_busco_list_${isolate_name}.tsv
		rm -r ${SAMPDATADIR}/run_${isolate_name}_BUSCO

		# Get end time of busco and calculate run time and append to time summary (and sum to total time used
		end=$SECONDS
		timeBUSCO=$((end - start))
		echo "BUSCO - ${timeBUSCO} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeBUSCO))
	# Prokka did not complete successfully and busco cant run (since it relies on prokka output)
	else
		echo "Prokka output not found, not able to process BUSCO"
	fi

	# Task: c-SSTAR for finding AR Genes
	write_Progress
	task_number=18
	echo "----- Running c-SSTAR for AR Gene identification -----"
	# c-SSTAR uses assembly and sample would have exited already if assembly did not complete, so no need to check
	# Get start time of ccstar
	start=$SECONDS

	# Run csstar in default mode from config.sh
	if [ ! -d "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/c-sstar${ResGANNCBI_srst2_filename}_${csstar_gapping}"
		mkdir -p "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}"
	fi
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES ${local_DBs}/singularities/cSSTAR.simg python3 /cSSTAR/c-SSTAR_${csstar_gapping}.py -g /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -s "${csim}" -d /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta --outdir /SAMPDIR/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping} > "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}.sstar"
	echo -e "c-SSTAR:1.1.01 -- python3 /cSSTAR/c-SSTAR_${csstar_gapping}.py -g ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta -s ${csim} -d ${local_DBs}/star/${ResGANNCBI_srst2_filename}_srst2.fasta --outdir ${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping} > ${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}.sstar\n" >> "${command_log_file}"
	###################################### FIND WAY TO CATCH FAILURE? !!!!!!!!!! ###############################

	# Goes through ResGANNCBI outfile and adds labels as well as resistance conferred to the beginning of the line
	# Takes .sstar file in and outputs as .sstar_grouped
	while IFS= read -r line; do

		#echo ${line}
		# Extract gene (label1) and allele (label2) from line, also force all characters to be lowercase
		label1=$(echo "${line}" | cut -d '	' -f3 | tr '[:upper:]' '[:lower:]')
		label2=$(echo "${line}" | cut -d '	' -f4 | tr '[:upper:]' '[:lower:]')
		# Determine what flags were thrown for this gene by csstar
		info1=""
		# Truncated allele
		if [[ "${label1}" = *"TRUNC" ]] && [[ "${label1}" != "str" ]]; then
			#echo "Label 1 was truncated"
			label1="${label1:0:${#label1} - 2}"
			info1="${info1}trunc-"
		fi
		# Likely novel allele
		if ( [[ "${label1}" = *"*"* ]] || [[ "${label1}" = *"*" ]] ) && [[ "${label1}" != "str" ]]; then
			#echo "Label 1 is likely novel"
			label1="${label1:0:${#label1} - 1}"
			info1="${info1}novel-"
		fi
		# Incomplete alignment length, Uncertainy exists in one allele
		if ( [[ "${label1}" = *"?"* ]] || [[ "${label1}" = *"?" ]] ) && [[ "${label1}" != "str" ]]; then
			#echo "Label 1 is uncertain due to incomplete alignment"
			label1="${label1:0:${#label1} - 1}"
			info1="${info1}alinc-"
		fi
		# Incomplete alignment length at edge
		if ( [[ "${label1}" = *"$"* ]] || [[ "${label1}" = *"$" ]] ) && [[ "${label1}" != "str" ]]; then
			#echo "Label 1 is uncertain due to incomplete alignment"
			label1="${label1:0:${#label1} - 1}"
			info1="${info1}edge-"
		fi
		# Removes character add-ons of genes and alleles, also lower cases all characters for searching later
		label1=$(echo "${label1,,}" | tr -d '*?$')
		label2=$(echo "${label2,,}" | tr -d '*?$')
		# Extract source database that AR gene match came from
		source=$(echo "${line,,}" | cut -d '	' -f1 | tr -d '[:space:]')
		# Extract the type of resistance that is conferred by the gene
		resistance=$(echo "${line}" | cut -d '	' -f2 | tr -d '[:space:]')
		# Trim contig identifier of spaces
		contig=$(echo "${line}" | cut -d '	' -f5 | tr -d '[:space:]')
		# Extract % from line
		percent=$(echo "${line}" | cut -d '	' -f6 | cut -d'%' -f1 | tr -d '[:space:]')
		# Determine length of query and subject sequences
		len1=$(echo "${line}" | cut -d '	' -f7 | tr -d '[:space:]')
		len2=$(echo "${line}" | cut -d '	' -f8 | tr -d '[:space:]')
		plen=$(echo "${line}" | cut -d '	' -f9 | tr -d '[:space:]')
		# Check and display any flags found, otherwise mark it as normal
		if [[ -z "${info1}" ]]; then
			info1="normal"
		else
			info1=${info1::-1}
		fi
		#printf "%-10s %-50s %-15s %-25s %-25s %-40s %-4s %-5d %-5d %-5d\\n" "${source}1" "${resistance}2" "${label1}3" "${info1}4" "${label2}5" "${contig}A" "${percent}B" "${len1}C" "${len2}D" "${plen}E"
		echo "${source}	${resistance}	${label1}	${info1}	${label2}	${contig}	${percent}	${len1}	${len2}	${plen}"
	done < "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}.sstar" > "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}.sstar_grouped"
	# Writes all AR genes to file based on %ID, %length, and finally length of gene
	sort -k7,7nr -k10,10nr -k8,8n "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}.sstar_grouped" > "${SAMPDATADIR}/c-sstar/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}_sstar_summary.txt"

	# Catches an empty or missing file, adding that no AMR genes were found if no file was created
	if [ ! -s "${SAMPDATADIR}/c-sstar/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}_sstar_summary.txt" ]; then
		echo "No anti-microbial genes were found using c-SSTAR with both resFinder and ARG-ANNOT DBs" > "${SAMPDATADIR}/c-sstar/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}_sstar_summary.txt"
	fi

	# Clean up
    # mv "${src}/${isolate_name}_scaffolds_trimmed"* "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}/"
	# mv "${src}/c-SSTAR_${isolate_name}_scaffolds_trimmed.log" "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}/"

	end=$SECONDS
	timestar=$((end - start))
	echo "c-SSTAR - ${timestar} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timestar))

	# Task: Run GAMMA on assembly
	write_Progress
	task_number=19
	echo "----- Running GAMMA -----"
	start=$SECONDS
	if [ ! -d "${SAMPDATADIR}/GAMMA" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/GAMMA"
		mkdir -p "${SAMPDATADIR}/GAMMA"
	fi
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/gamma:1.4--hdfd78af_0 GAMMA.py /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta /SAMPDIR/GAMMA/${isolate_name}.${ResGANNCBI_srst2_filename}
	echo -e "GAMMA:1.4 -- GAMMA.py ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta ${local_DBs}/star/${ResGANNCBI_srst2_filename}_srst2.fasta ${SAMPDATADIR}/GAMMA/${isolate_name}.${ResGANNCBI_srst2_filename}" >> "${command_log_file}"

	python3 ${src}/GAMMA_ResGANNCBI_file_converter.py ${SAMPDATADIR}/GAMMA/${isolate_name}.${ResGANNCBI_srst2_filename}.gamma

	end=$SECONDS
	timeGAMMA=$((end - start))
	echo "GAMMA - ${timeGAMMA} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeGAMMA))

	# Task: Get MLST profile
	write_Progress
	task_number=20
	echo "----- Running MLST -----"
	if [ ! -d "${SAMPDATADIR}/MLST" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/MLST"
		mkdir -p "${SAMPDATADIR}/MLST"
	fi
	start=$SECONDS
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/mlst:2.19.0--hdfd78af_1 mlst /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta > "${SAMPDATADIR}/MLST/${isolate_name}.mlst"
	echo -e "mlst:2.19 -- mlst ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta > ${SAMPDATADIR}/MLST/${isolate_name}.mlst\n" >> "${command_log_file}"
	python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}.mlst" -t standard -d ${local_DBs}/pubmlst
	type=$(head -n1 ${SAMPDATADIR}/MLST/${isolate_name}.mlst | cut -d' ' -f3)
	if [[ "${genus}_${species}" = "Acinetobacter_baumannii" ]]; then
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/mlst:2.19.0--hdfd78af_1 mlst --scheme "abaumannii" /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta > "${SAMPDATADIR}/MLST/${isolate_name}_abaumannii.mlst"
		echo -e "mlst:2.19 -- mlst --scheme \"abaumannii\" ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta > ${SAMPDATADIR}/MLST/${isolate_name}_abaumannii.mlst\n" >> "${command_log_file}"
		python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_abaumannii.mlst" -t standard -d ${local_DBs}/pubmlst
		mv "${SAMPDATADIR}/MLST/${isolate_name}_abaumannii.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Oxford.mlst"
		mv "${SAMPDATADIR}/MLST/${isolate_name}.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Pasteur.mlst"
		#Check for "-", unidentified type
		type1=$(tail -n1 ${SAMPDATADIR}/MLST/${isolate_name}_Oxford.mlst | cut -d' ' -f3)
		type2=$(head -n1 ${SAMPDATADIR}/MLST/${isolate_name}_Pasteur.mlst | cut -d' ' -f3)
		if [[ "${type1}" = "-" ]]; then
			singularity -s exec ${local_DBs}/singularities/srst2.simg python /usr/local/bin/getmlst.py --species "Acinetobacter baumannii#1" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
			echo -e "srst2:0.2.0 -- getmlst.py --species \"Acinetobacter baumannii#1\" > ${SAMPDATADIR}/MLST/srst2/getmlst.out\n" >> "${command_log_file}"
			sed -i -e 's/Oxf_//g' "${SAMPDATADIR}/MLST/srst2/Acinetobacter_baumannii#1.fasta"
			sed -i -e 's/Oxf_//g' "${SAMPDATADIR}/MLST/srst2/abaumannii.txt"
			db_name="Oxford"
			suggested_command=$(tail -n2 "${SAMPDATADIR}/MLST/srst2/getmlst.out" | head -n1)
			mlst_db=$(echo "${suggested_command}" | cut -d' ' -f11)
			mlst_defs=$(echo "${suggested_command}" | cut -d' ' -f13)
			mlst_delimiter=$(echo "${suggested_command}" | cut -d' ' -f15)
			if [[ "${mlst_delimiter}" != "'_'" ]]; then
				echo "Unknown delimiter - \"${mlst_delimiter}\""
			else
				mlst_delimiter="_"
			fi
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${local_DBs}/singularities/srst2.simg python /usr/local/bin/srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}
			echo -e "srst2:0.2.0 -- srst2 --input_pe ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output ${SAMPDATADIR}/srst2/MLST/srst2/${isolate_name} --mlst_db ${SAMPDATADIR}/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}\n" >> "${command_log_file}"
			today=$(date "+%Y-%m-%d")

			# Cleans up extra files and renames output file
			mv "${SAMPDATADIR}/MLST/srst2/${isolate_name}__mlst__Acinetobacter_baumannii#1__results.txt" "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Acinetobacter_baumannii#1-${db_name}.mlst"
			mv "${SAMPDATADIR}/MLST/srst2/mlst_data_download_Acinetobacter_baumannii#1_${today}.log" "${SAMPDATADIR}/MLST/"
			rm -r "${SAMPDATADIR}/MLST/srst2"

			if [[ -f "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.Acinetobacter_baumannii#1.pileup" ]]; then
				rm -r "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.Acinetobacter_baumannii#1.pileup"
			fi
			if [[ -f "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.Acinetobacter_baumannii#1.sorted.bam" ]]; then
				rm -r "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.Acinetobacter_baumannii#1.sorted.bam"
			fi
			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Acinetobacter_baumannii#1-${db_name}.mlst" -t srst2 -d ${local_DBs}/pubmlst
		fi
		if [[ "${type2}" = "-" ]]; then
			singularity -s exec ${local_DBs}/singularities/srst2.simg python /usr/local/bin/getmlst.py --species "Acinetobacter baumannii#2" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
			echo -e "srst2:0.2.0 -- getmlst.py --species \"Acinetobacter baumannii#2\" > ${SAMPDATADIR}/MLST/srst2/getmlst.out\n" >> "${command_log_file}"
			sed -i -e 's/Pas_//g' "${SAMPDATADIR}/MLST/srst2/Acinetobacter_baumannii#2.fasta"
			sed -i -e 's/Pas_//g' "${SAMPDATADIR}/MLST/srst2/abaumannii_2.txt"
			db_name="Pasteur"
			suggested_command=$(tail -n2 "${SAMPDATADIR}/MLST/srst2/getmlst.out" | head -n1)
			mlst_db=$(echo "${suggested_command}" | cut -d' ' -f11)
			mlst_defs=$(echo "${suggested_command}" | cut -d' ' -f13)
			mlst_delimiter=$(echo "${suggested_command}" | cut -d' ' -f15)
			if [[ "${mlst_delimiter}" != "'_'" ]]; then
				echo "Unknown delimiter - \"${mlst_delimiter}\""
			else
				mlst_delimiter="_"
			fi
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${local_DBs}/singularities/srst2.simg python /usr/local/bin/srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}
			echo -e "srst2:0.2.0 -- srst2 --input_pe ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output ${SAMPDATADIR}/srst2/MLST/srst2/${isolate_name} --mlst_db ${SAMPDATADIR}/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}\n" >> "${command_log_file}"
			today=$(date "+%Y-%m-%d")

			# Cleans up extra files and renames output file
			mv "${SAMPDATADIR}/MLST/srst2/${isolate_name}__mlst__Acinetobacter_baumannii#2__results.txt" "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Acinetobacter_baumannii#2-${db_name}.mlst"
			mv "${SAMPDATADIR}/MLST/srst2/mlst_data_download_Acinetobacter_baumannii#2_${today}.log" "${SAMPDATADIR}/MLST/"
			rm -r "${SAMPDATADIR}/MLST/srst2"

			if [[ -f "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.Acinetobacter_baumannii#2.pileup" ]]; then
				rm -r "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.Acinetobacter_baumannii#2.pileup"
			fi
			if [[ -f "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.Acinetobacter_baumannii#2.sorted.bam" ]]; then
				rm -r "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.Acinetobacter_baumannii#2.sorted.bam"
			fi
			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Acinetobacter_baumannii#2-${db_name}.mlst" -t srst2 -d ${local_DBs}/pubmlst
		fi
	elif [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
		# Verify that ecoli_2 is default and change accordingly
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/mlst:2.19.0--hdfd78af_1 mlst --scheme "ecoli_2" /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta > "${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst"
		echo -e "mlst:2.19 -- mlst --scheme \"ecoli_2\" ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta > ${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst\n" >> "${command_log_file}"
		python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst" -t standard -d ${local_DBs}/pubmlst
		mv "${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Pasteur.mlst"
		mv "${SAMPDATADIR}/MLST/${isolate_name}.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Achtman.mlst"
		type2=$(tail -n1 ${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst | cut -d' ' -f3)
		type1=$(head -n1 ${SAMPDATADIR}/MLST/${isolate_name}.mlst | cut -d' ' -f3)
		if [[ "${type1}" = "-" ]]; then
			singularity -s exec ${local_DBs}/singularities/srst2.simg python /usr/local/bin/getmlst.py --species "Escherichia coli#1" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
			echo -e "srst2:0.2.0 -- getmlst.py --species \"Escherichia coli#1\" > ${SAMPDATADIR}/MLST/srst2/getmlst.out\n" >> "${command_log_file}"
			db_name="Achtman"
			suggested_command=$(tail -n2 "${SAMPDATADIR}/MLST/srst2/getmlst.out" | head -n1)
			mlst_db=$(echo "${suggested_command}" | cut -d' ' -f11)
			mlst_defs=$(echo "${suggested_command}" | cut -d' ' -f13)
			mlst_delimiter=$(echo "${suggested_command}" | cut -d' ' -f15)
			if [[ "${mlst_delimiter}" != "'_'" ]]; then
				echo "Unknown delimiter - \"${mlst_delimiter}\""
			else
				mlst_delimiter="_"
				#echo "Delimiter is OK (${mlst_delimiter})"
			fi
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${local_DBs}/singularities/srst2.simg python /usr/local/bin/srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}
			echo -e "srst2:0.2.0 -- srst2 --input_pe ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output ${SAMPDATADIR}/srst2/MLST/srst2/${isolate_name} --mlst_db ${SAMPDATADIR}/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}\n" >> "${command_log_file}"
			today=$(date "+%Y-%m-%d")
					# Cleans up extra files and renames output file
			mv "${SAMPDATADIR}/MLST/srst2/${isolate_name}__mlst__Escherichia_coli#1__results.txt" "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Escherichia_coli#1-${db_name}.mlst"
			mv "${SAMPDATADIR}/MLST/srst2/mlst_data_download_Escherichia_coli#1_${today}.log" "${SAMPDATADIR}/MLST/"
			rm -r "${SAMPDATADIR}/MLST/srst2"

			if [[ -f "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.Escherichia_coli#1.pileup" ]]; then
				rm -r "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.Escherichia_coli#1.pileup"
			fi
			if [[ -f "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.Escherichia_coli#1.sorted.bam" ]]; then
				rm -r "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.Escherichia_coli#1.sorted.bam"
			fi

			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Escherichia_coli#1-${db_name}.mlst" -t srst2 -d ${local_DBs}/pubmlst
		fi
		if [[ "${type2}" = "-" ]]; then
			singularity -s exec ${local_DBs}/singularities/srst2.simg python /usr/local/bin/getmlst.py --species "Escherichia coli#2" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
			echo -e "srst2:0.2.0 -- getmlst.py --species \"Escherichia coli#2\" > ${SAMPDATADIR}/MLST/srst2/getmlst.out\n" >> "${command_log_file}"
			db_name="Pasteur"
			suggested_command=$(tail -n2 "${SAMPDATADIR}/MLST/srst2/getmlst.out" | head -n1)
			mlst_db=$(echo "${suggested_command}" | cut -d' ' -f11)
			mlst_defs=$(echo "${suggested_command}" | cut -d' ' -f13)
			mlst_delimiter=$(echo "${suggested_command}" | cut -d' ' -f15)
			if [[ "${mlst_delimiter}" != "'_'" ]]; then
				echo "Unknown delimiter - \"${mlst_delimiter}\""
			else
				mlst_delimiter="_"
				#echo "Delimiter is OK (${mlst_delimiter})"
			fi
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${local_DBs}/singularities/srst2.simg python /usr/local/bin/srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}
			echo -e "srst2:0.2.0 -- srst2 --input_pe ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output ${SAMPDATADIR}/srst2/MLST/srst2/${isolate_name} --mlst_db ${SAMPDATADIR}/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}\n" >> "${command_log_file}"
			today=$(date "+%Y-%m-%d")
					# Cleans up extra files and renames output file
			mv "${SAMPDATADIR}/MLST/srst2/${isolate_name}__mlst__Escherichia_coli#2__results.txt" "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Escherichia_coli#2-${db_name}.mlst"
			mv "${SAMPDATADIR}/MLST/srst2/mlst_data_download_Escherichia_coli#2_${today}.log" "${SAMPDATADIR}/MLST/"
			rm -r "${SAMPDATADIR}/MLST/srst2"

			if [[ -f "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.Escherichia_coli#2.pileup" ]]; then
				rm -r "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.Escherichia_coli#2.pileup"
			fi
			if [[ -f "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.Escherichia_coli#2.sorted.bam" ]]; then
				rm -r "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.Escherichia_coli#2.sorted.bam"
			fi

			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Escherichia_coli#2-${db_name}.mlst" -t srst2 -d ${local_DBs}/pubmlst
		fi
	else
		if [[ "${type}" == "-" ]]; then
			singularity -s exec ${local_DBs}/singularities/srst2.simg python /usr/local/bin/getmlst.py --species "Escherichia coli#1" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
			echo -e "srst2:0.2.0 -- getmlst.py --species \"Escherichia coli#1\" > ${SAMPDATADIR}/MLST/srst2/getmlst.out\n" >> "${command_log_file}"
			db_name="Pasteur"
			suggested_command=$(tail -n2 "${SAMPDATADIR}/MLST/srst2/getmlst.out" | head -n1)
			mlst_db=$(echo "${suggested_command}" | cut -d' ' -f11)
			mlst_defs=$(echo "${suggested_command}" | cut -d' ' -f13)
			mlst_delimiter=$(echo "${suggested_command}" | cut -d' ' -f15)
			if [[ "${mlst_delimiter}" != "'_'" ]]; then
				echo "Unknown delimiter - \"${mlst_delimiter}\""
			else
				mlst_delimiter="_"
				#echo "Delimiter is OK (${mlst_delimiter})"
			fi
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${local_DBs}/singularities/srst2.simg python /usr/local/bin/srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}
			echo -e "srst2:0.2.0 -- srst2 --input_pe ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output ${SAMPDATADIR}/srst2/MLST/srst2/${isolate_name} --mlst_db ${SAMPDATADIR}/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}\n" >> "${command_log_file}"
			today=$(date "+%Y-%m-%d")
					# Cleans up extra files and renames output file
			mv "${SAMPDATADIR}/MLST/srst2/${isolate_name}__mlst__${genus}_${species}__results.txt" "${SAMPDATADIR}/MLST/${isolate_name}_srst2_${genus}_${species}-${db_name}.mlst"
			mv "${SAMPDATADIR}/MLST/srst2/mlst_data_download_${genus}_${species}_${today}.log" "${SAMPDATADIR}/MLST/"
			rm -r "${SAMPDATADIR}/MLST/srst2"

			if [[ -f "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.${genus}_${species}.pileup" ]]; then
				rm -r "${SAMPDATADIR}/MLST/srst2/${isolate_name}__${isolate_name}.${genus}_${species}.pileup"
			fi
			if [[ -f "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.${genus}_${species}.sorted.bam" ]]; then
				rm -r "${SAMPDATADIR}/MLST/${isolate_name}__${isolate_name}.${genus}_${species}.sorted.bam"
			fi
			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_${genus}_${species}-${db_name}.mlst" -t srst2 -d ${local_DBs}/pubmlst
		fi
		mv "${SAMPDATADIR}/MLST/${isolate_name}.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Pasteur.mlst"
	fi
	end=$SECONDS
	timeMLST=$((end - start))
	echo "MLST - ${timeMLST} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeMLST))

	# Task: Try to find any plasmids
	write_Progress
	task_number=21
	echo "----- Identifying plasmid replicons using plasmidFinder -----"
	start=$SECONDS
	if [[ ! -d "${SAMPDATADIR}/plasmidFinder" ]]; then
		mkdir "${SAMPDATADIR}/plasmidFinder"
	fi
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${local_DBs}/singularities/plasmidFinder_with_DB.simg plasmidfinder.py -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmidFinder -p /opt/plasmidfinder_db -t ${plasmidFinder_identity}
	echo -e "plasmidFinder:2.1 -- plasmidfinder.py -i ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o ${SAMPDATADIR}/plasmidFinder -p /opt/plasmidfinder_db -t ${plasmidFinder_identity}\n" >> "${command_log_file}"
	python "${src}/json_plasmidFinder_converter.py" -i "${SAMPDATADIR}/plasmidFinder/data.json" -o "${SAMPDATADIR}/plasmidFinder/${isolate_name}_results_table_summary.txt"

	end=$SECONDS
	timeplasfin=$((end - start))
	echo "plasmidFinder - ${timeplasfin} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeplasfin))

	# # Run plasFlow if isolate is from the Enterobacteriaceae family  ##### When should we check if this will be expanded?
	# if [[ "${family}" == "Enterobacteriaceae" ]]; then
	# 	start=$SECONDS
	# 	# Create output directory
	# 	if [[ ! -d "${SAMPDATADIR}/plasFlow" ]]; then
	# 		mkdir "${SAMPDATADIR}/plasFlow"
	# 	fi
	# 	# Task: Run all tools associated with reconsructing plasmids using plasFlow
	# 	write_Progress
	# 	task_number=22
	# 	if [[ -s "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" ]]; then
	# 		# Trim contigs a little to 2000 and larger and put through plasflow.
	# 		python3 "${src}/removeShortContigs.py" -i "${SAMPDATADIR}/Assembly/scaffolds.fasta" -t 2000 -s "normal_SPAdes"
	# 		# Renames headers of fasta files accordingly
	# 		python3 "${src}/fasta_headers.py" -i "${SAMPDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta" -o "${SAMPDATADIR}/plasFlow/${isolate_name}_scaffolds_trimmed_2000.fasta"
	# 		# Removes intermeidate fasta file
	# 		rm -r "${SAMPDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta"
	# 		# Run plasflow on newly trimmed assembly file
	# 		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasflow:1.1.0--py35_0 PlasFlow.py --input /SAMPDIR/plasFlow/${isolate_name}_scaffolds_trimmed_2000.fasta --output /SAMPDIR/plasFlow/${isolate_name}_plasFlow_results.tsv --threshold 0.7
	# 		echo -e "PlasFlow:1.1.0 -- PlasFlow.py --input ${SAMPDATADIR}/plasFlow/${isolate_name}_scaffolds_trimmed_2000.fasta --output ${SAMPDATADIR}/plasFlow/${isolate_name}_plasFlow_results.tsv --threshold 0.7\n" >> "${command_log_file}"
	# 		mkdir ${SAMPDATADIR}/plasFlow/bowtie2-index/
	# 		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/bowtie2:2.4.4--py37h13ad519_0 bowtie2-build -f /SAMPDIR/plasFlow/${isolate_name}_plasFlow_results.tsv_chromosomes.fasta /SAMPDIR/plasFlow/bowtie2-index/bowtie2_${isolate_name}_chr
	# 		echo -e "bowtie2:2.4.4 -- bowtie2-build -f ${SAMPDATADIR}/plasFlow/${isolate_name}_plasFlow_results.tsv_chromosomes.fasta ${SAMPDATADIR}/plasFlow/bowtie2-index/bowtie2_${isolate_name}_chr\n" >> "${command_log_file}"
	# 		mkdir ${SAMPDATADIR}/plasFlow/filtered_reads_70/
	# 		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/bowtie2:2.4.4--py37h13ad519_0 bowtie2 -x /SAMPDIR/plasFlow/bowtie2-index/bowtie2_${isolate_name}_chr -1 /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq -2 /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq -S /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}.sam -p ${procs} --local
	# 		echo -e "bowtie2:2.4.4 -- bowtie2 -x ${SAMPDATADIR}/plasFlow/bowtie2-index/bowtie2_${isolate_name}_chr -1 ${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq -2 ${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq -S ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}.sam -p ${procs} --local\n" >> "${command_log_file}"
	# 		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/samtools:1.14--hb421002_0 samtools view -bS /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}.sam > "${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}.bam"
	# 		echo -e "samtools:1.14 -- samtools view -bS ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}.sam > ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}.bam\n" >> "${command_log_file}"
	# 		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/samtools:1.14--hb421002_0 samtools sort -n /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}.bam -o /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}.bam.sorted
	# 		echo -e "samtools:1.14 -- samtools sort -n ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}.bam -o ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}.bam.sorted\n" >> "${command_log_file}"
	# 		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/bedtools:2.30.0--h7d7f7ad_2 bamToFastq -i /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}.bam.sorted -fq /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}_R1_bacterial.fastq -fq2 /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}_R2_bacterial.fastq
	# 		echo -e "bedtools:2.30.0 -- bamToFastq -i ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}.bam.sorted -fq ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}_R1_bacterial.fastq -fq2 ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}_R2_bacterial.fastq\n" >> "${command_log_file}"
	# 		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/unicycler:0.4.4--py37h8b12597_2 unicycler -1 /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}_R1_bacterial.fastq -2 /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}_R2_bacterial.fastq -o /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly -t "${procs}"
	# 		echo -e "unicycler:0.4.4 -- unicycler -1 ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}_R1_bacterial.fastq -2 ${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}_R2_bacterial.fastq -o ${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly -t ${procs}\n" >> "${command_log_file}"
	# 		mv "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/assembly.fasta" "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_original.fasta"
	# 		mv "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/assembly.gfa" "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_assembly.gfa"
	# 		python3 "${src}/fasta_headers_plasFlow.py" -i "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_original.fasta" -o "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly.fasta"
	# 		python3 "${src}/removeShortContigs.py" -i "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly.fasta" -t 500 -s "plasFlow"
	# 		mv "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly.fasta.TRIMMED.fasta" "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta"
	# 		rm "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly.fasta"
	# 	else
	# 		echo "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta (Assembly) not found, cant do anything"
	# 	fi
	#
	# 	# Check quality of plasFlow/Unicycler assembly
	# 	write_Progress
	# 	task_number=23
	# 	if [[ -s "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta" ]]; then
	# 		if [ ! -d "${SAMPDATADIR}/Assembly_Stats_plasFlow" ]; then
	# 			echo "Creating ${SAMPDATADIR}/Assembly_Stats_plasFlow"
	# 	 		mkdir -p "${SAMPDATADIR}/Assembly_Stats_plasFlow"
	# 	 	fi
	# 	 	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/quast:5.0.2--py37pl5262hfecc14a_5 quast.py --no-icarus --no-html --no-snps -o /SAMPDIR/Assembly_Stats_plasFlow /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta
	# 		echo -e "QUAST:5.0.2-- python3 quast.py -o ${SAMPDATADIR}/Assembly_Stats_plasFlow ${SAMPDATADIR}/Assembly/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta\n" >> "${command_log_file}"
	# 		mv "${SAMPDATADIR}/Assembly_Stats_plasFlow/report.txt" "${SAMPDATADIR}/Assembly_Stats_plasFlow/${isolate_name}_report.txt"
	# 	 	mv "${SAMPDATADIR}/Assembly_Stats_plasFlow/report.tsv" "${SAMPDATADIR}/Assembly_Stats_plasFlow/${isolate_name}_report.tsv"
	# 	 else
	# 		echo "No plasFlow assembly (${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta)"
	# 	fi
	#
	# 	# Task: Run c-SSTAR on plasFlow Assembly
	# 	write_Progress
	# 	task_number=24
	# 	# Creates the output c-sstar folder if it does not exist yet
	# 	if [ ! -d "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}" ]; then  #create outdir if absent
	# 		echo "Creating ${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}"
	# 		mkdir -p "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}"
	# 	fi
	# 	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES ${local_DBs}/singularities/cSSTAR.simg python3 /cSSTAR/c-SSTAR_${csstar_gapping}.py -g /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta -s "${cpsim}" -d /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta -o /SAMPDIR/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping} > "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar"
	# 	echo -e "c-SSTAR:1.1.01 -- python3 /cSSTAR/c-SSTAR_${csstar_gapping}.py -g ${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta -s ${cpsim} -d /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta > ${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar\n" >> "${command_log_file}"
	#
	#
	# 	###################################### FIND WAY TO CATCH FAILURE !!!!!!!!!! ###############################
	#
	# 	# Goes through ResGANNCBI outfile and adds labels as well as resistance conferred to the beginning of the line
	# 	# Takes .sstar file in and outputs as .sstar_grouped
	# 	while IFS= read -r line; do
	#
	# 		#echo ${line}
	# 		# Extract gene (label1) and allele (label2) from line, also force all characters to be lowercase
	# 		label1=$(echo "${line}" | cut -d '	' -f3 | tr '[:upper:]' '[:lower:]')
	# 		label2=$(echo "${line}" | cut -d '	' -f4 | tr '[:upper:]' '[:lower:]')
	# 		# Determine what flags were thrown for this gene by csstar
	# 		info1=""
	# 		# Truncated allele
	# 		if [[ "${label1}" = *"TRUNC" ]] && [[ "${label1}" != "str" ]]; then
	# 			#echo "Label 1 was truncated"
	# 			label1="${label1:0:${#label1} - 2}"
	# 			info1="${info1}trunc-"
	# 		fi
	# 		# Likely novel allele
	# 		if ( [[ "${label1}" = *"*"* ]] || [[ "${label1}" = *"*" ]] ) && [[ "${label1}" != "str" ]]; then
	# 			#echo "Label 1 is likely novel"
	# 			label1="${label1:0:${#label1} - 1}"
	# 			info1="${info1}novel-"
	# 		fi
	# 		# Incomplete alignment length, Uncertainy exists in one allele
	# 		if ( [[ "${label1}" = *"?"* ]] || [[ "${label1}" = *"?" ]] ) && [[ "${label1}" != "str" ]]; then
	# 			#echo "Label 1 is uncertain due to incomplete alignment"
	# 			label1="${label1:0:${#label1} - 1}"
	# 			info1="${info1}alinc-"
	# 		fi
	# 		# Incomplete alignment length at edge
	# 		if ( [[ "${label1}" = *"$"* ]] || [[ "${label1}" = *"$" ]] ) && [[ "${label1}" != "str" ]]; then
	# 			#echo "Label 1 is uncertain due to incomplete alignment"
	# 			label1="${label1:0:${#label1} - 1}"
	# 			info1="${info1}edge-"
	# 		fi
	# 		# Removes character add-ons of genes and alleles, also lower cases all characters for searching later
	# 		label1=$(echo "${label1,,}" | tr -d '*?$')
	# 		label2=$(echo "${label2,,}" | tr -d '*?$')
	# 		# Extract source database that AR gene match came from
	# 		source=$(echo "${line,,}" | cut -d '	' -f1 | tr -d '[:space:]')
	# 		# Extract the type of resistance that is conferred by the gene
	# 		resistance=$(echo "${line}" | cut -d '	' -f2 | tr -d '[:space:]')
	# 		# Trim contig identifier of spaces
	# 		contig=$(echo "${line}" | cut -d '	' -f5 | tr -d '[:space:]')
	# 		# Extract % from line
	# 		percent=$(echo "${line}" | cut -d '	' -f6 | cut -d'%' -f1 | tr -d '[:space:]')
	# 		# Determine length of query and subject sequences
	# 		len1=$(echo "${line}" | cut -d '	' -f7 | tr -d '[:space:]')
	# 		len2=$(echo "${line}" | cut -d '	' -f8 | tr -d '[:space:]')
	# 		plen=$(echo "${line}" | cut -d '	' -f9 | tr -d '[:space:]')
	#
	# 		# Check and display any flags found, otherwise mark it as normal
	# 		if [[ -z "${info1}" ]]; then
	# 			info1="normal"
	# 		else
	# 			info1=${info1::-1}
	# 		fi
	# 		#printf "%-10s %-50s %-15s %-25s %-25s %-40s %-4s %-5d %-5d %-5d\\n" "${source}1" "${resistance}2" "${label1}3" "${info1}4" "${label2}5" "${contig}A" "${percent}B" "${len1}C" "${len2}D" "${plen}E"
	# 		echo "${source}	${resistance}	${label1}	${info1}	${label2}	${contig}	${percent}	${len1}	${len2}	${plen}"
	# 	done < "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar" > "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar_grouped"
	# 	# Writes all AR genes to file based on %ID, %length, and finally length of gene
	# 	sort -k7,7nr -k10,10nr -k8,8n "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar_grouped" > "${SAMPDATADIR}/c-sstar_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}_sstar_summary.txt"
	#
	# 	# Catches an empty or missing file, adding that no AMR genes were found if no file was created
	# 	if [ ! -s "${SAMPDATADIR}/c-sstar_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}_sstar_summary.txt" ]; then
	# 		echo "No anti-microbial genes were found using c-SSTAR with both resFinder and ARG-ANNOT DBs" > "${SAMPDATADIR}/c-sstar_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}_sstar_summary.txt"
	# 	fi
	#
	# 	# Clean up
	# 	mv "${src}/${isolate_name}_plasmid_assembly_trimmed"* "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/"
	# 	mv "${src}/c-SSTAR_${isolate_name}_plasmid_assembly_trimmed.log" "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/"
	#
	#
	# 	# Task: Try to find any plasmid replicons on plasFlow assembly
	# 	write_Progress
	# 	task_number=25
	# 	echo "----- Identifying plasmids using plasmidFinder -----"
	# 	if [[ ! -d "${SAMPDATADIR}/plasmidFinder_on_plasFlow" ]]; then
	# 		mkdir "${SAMPDATADIR}/plasmidFinder_on_plasFlow"
	# 	fi
	# 	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${local_DBs}/singularities/plasmidFinder_with_DB.simg plasmidfinder.py -i /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta -o /SAMPDIR/plasmidFinder_on_plasFlow -p /opt/plasmidfinder_db -t ${plasmidFinder_identity}
	# 	echo -e "plasmidFinder:2.1 -- plasmidfinder.py -i ${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta -o ${SAMPDATADIR}/plasmidFinder_on_plasFlow -p /opt/plasmidfinder_db -t ${plasmidFinder_identity}\n" >> "${command_log_file}"
	# 	python "${src}/json_plasmidFinder_converter.py" -i "${SAMPDATADIR}/plasmidFinder_on_plasFlow/data.json" -o "${SAMPDATADIR}/plasmidFinder_on_plasFlow/${isolate_name}_results_table_summary.txt"
	#
	# 	# Task: Use GAMMA to find any AR genes in plasFlow assembly_length
	# 	write_Progress
	# 	task_number=26
	# 	echo "----- Identifying AR genes with GAMMA -----"
	# 	if [ ! -d "${SAMPDATADIR}/GAMMA_plasFlow" ]; then  #create outdir if absent
	# 		echo "Creating ${SAMPDATADIR}/GAMMA_plasFlow"
	# 		mkdir -p "${SAMPDATADIR}/GAMMA_plasFlow"
	# 	fi
	# 	#singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES ${local_DBs}/singularities/GAMMA_quaisar.simg python3 /GAMMA/GAMMA_quaisar.py -i /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta -d /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta -o /SAMPDIR/GAMMA_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.GAMMA
	# 	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES	docker://quay.io/biocontainers/gamma:1.4--hdfd78af_0 GAMMA.py /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta /SAMPDIR/GAMMA_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}
	#
	# 	#echo -e "GAMMA:4.7.4 -- python3 /GAMMA/GAMMA_quaisar.py -i ${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta -d /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta -o ${SAMPDATADIR}/GAMMA_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.GAMMA\n" >> "${command_log_file}"
	# 	echo -e "GAMMA:1.4 -- GAMMA.py ${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta ${SAMPDATADIR}/GAMMA_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}\n" >> "${command_log_file}"
	#
	# 	python3 ${src}/GAMMA_ResGANNCBI_file_converter.py ${SAMPDIR}/GAMMA_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.gamma
	#
	# 	end=$SECONDS
	# 	timeplasflow=$((end - start))
	# 	echo "plasFlow - ${timeplasflow} seconds" >> "${time_summary}"
	# 	totaltime=$((totaltime + timeplasflow))
	# fi

	# Task: Create stats file
	write_Progress
	task_number=27
	"${src}/validate_piperun.sh" "${SAMPDATADIR}" "${local_DBs}" "${src}" "${csstar_gapping}" "${csim}" "${cpsim}" > "${SAMPDATADIR}/${isolate_name}_pipeline_stats.txt"

	# Task: Clean sample folder
	write_Progress
	task_number=28
	"${src}/sample_cleaner.sh" "${SAMPDATADIR}"

	rm -r "${src}/tmp"

	# Extra dump cleanse in case anything else failed
	if [ -n "$(find "${src}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
		echo "Found core dump files at end of processing ${isolate_name} and attempting to delete"
		find "${src}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
	fi

	global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

	# Append total time to bottom of time summary
	echo "Total time: ${totaltime} seconds" >> "${time_summary}"
	echo "Completed at ${global_end_time}"

	echo -e "\n${isolate} complete at ${global_end_time}\n" >> "${log_file}"

	# Designate end of this sample #
	echo "

					End of sample ${isolate_name}
					completed at ${global_end_time}

	"
	echo -e "\n\n" >> "${command_log_file}"
	loop_inc=$(( loop_inc + 1 ))
done

# Task: Concatenates lists if this run was an addition to an already processed folder
write_Progress
run_task_id=7
if [[ -f "${PROJDATADIR}/${PROJECT}_list_original.txt" ]]; then
	cat "${PROJDATADIR}/${PROJECT}_list.txt" >> "${PROJDATADIR}/${PROJECT}_list_original.txt"
	rm "${PROJDATADIR}/${PROJECT}_list.txt"
	mv "${PROJDATADIR}/${PROJECT}_list_original.txt" "${PROJDATADIR}/${PROJECT}_list.txt"
fi

# Task: Run the Seqlog creator on the proper file
#write_Progress
#run_task_id=8
# declare -A mmb_bugs
# while IFS= read -r bug_lines  || [ -n "$bug_lines" ]; do
# 	bug_genus=$(echo "${bug_lines}" | cut -d'	' -f1)
# 	bug_species=$(echo "${bug_lines}" | cut -d'	' -f2)
# 	bug_info=$(echo "${bug_lines}" | cut -d'	' -f3-)
# 	bug_size=$(echo "${bug_lines}" | cut -d'	' -f6)
# 	bug_name="${bug_genus:0:1}.${bug_species}"
# 	#echo "Should be adding ${bug_size} for ${bug_name}"
# 	mmb_bugs["${bug_name}"]="${bug_size}"
# done < ${local_DBs}/MMB_Bugs.txt

# Set output folder as directory of input list
> "${PROJDATADIR}/Seqlog_output.txt"

# Task: Create Seqlog info
write_Progress
run_task_id=9

# Create header for file
echo "Isolate_ID	Analysis_Date	KRAKEN ID Raw Reads	KRAKEN ID Assembly	16S BLAST ID	Q20_Total_[bp]	Q30_Total_[bp]	Q20_R1_[bp]	Q20_R2_[bp]	Q20_R1_[%]	Q20_R2_[%]	Q30_R1_[bp]	Q30_R2_[bp]	Q30_R1_[%]	Q30_R2_[%]	Total_Sequenced_[bp]	Total_Sequenced_[reads]	Estimated coverage	Contigs	Cumulative_Length_Assembly (bp)	Assembly_Ratio	BUSCO	ANI" > "${PROJDATADIR}/Seqlog_output.txt"
while IFS= read -r var || [ -n "$var" ]; do
	# Current (12/17/18) order of expected run output
	#  kraken - QC - estimated coverage - #contigs - cumulative length assembly - BUSCO - ANI
	isolate_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	SAMPDATADIR="${PROJDATADIR}/${isolate_name}"

	# Creates default values in case they are not filled in later
	g_s_assembled="Unidentified"
	genus_post="not_assigned"
	species_post="not_assigned"
	# Pulls species and genus_post information from kraken out of assembly
	if [[ -s "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_kraken_summary_assembled_BP.txt" ]]; then
		while IFS= read -r line; do
			first=${line::1}
			if [ "${first}" = "s" ]
			then
				species_post=$(echo "${line}" | awk -F ' ' '{print $4}')
			elif [ "${first}" = "G" ]
			then
				genus_post=$(echo "${line}" | awk -F ' ' '{print $4}')
			fi
		done < "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_kraken_summary_assembled_BP.txt"
		g_s_assembled="${genus_post} ${species_post}"
		#echo "${g_s_assembly}"
	elif [[ ! -f "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" ]]; then
		#echo "Cant find ${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta"
		g_s_assembled="Failed_Assembly"
	fi

	g_s_reads="Unidentified"
	genus_reads="not_assigned"
	species_reads="not_assigned"
	# Pulls species and genus_16s information from kraken out of reads
	if [[ -s "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_kraken_summary_paired.txt" ]]; then
		while IFS= read -r line; do
			first=${line::1}
			if [ "${first}" = "s" ]
			then
				species_reads=$(echo "${line}" | awk -F ' ' '{print $4}')
			elif [ "${first}" = "G" ]
			then
				genus_reads=$(echo "${line}" | awk -F ' ' '{print $4}')
			fi
		done < "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_kraken_summary_paired.txt"
		g_s_reads="${genus_reads} ${species_reads}"
		#echo "${g_s}"
	fi
	# Pulls species and genus_16s information from 16s
	g_s_16s="Unidentified"
	genus_16s="not_assigned"
	species_16s="not_assigned"
	if [[ -s "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt" ]]; then
		info=$(tail -n 1 "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt")
		type=$(echo "${info}" | cut -d' ' -f1)
		genus_16s=$(echo "${info}" | cut -d'	' -f3 | cut -d' ' -f1)
		species_16s=$(echo "${info}" | cut -d'	' -f3 | cut -d' ' -f2)
		if [[ "${genus_16s}" = "No_16s_sequences_found" ]] && [[ "${genus_16s}" = "No_16s_sequences_found" ]]; then
			g_s_16s="${genus_16s}"
		else
			g_s_16s="${genus_16s} ${species_16s}"
		fi
	elif [[ ! -f "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" ]]; then
		g_s_16s="Failed_Assembly"
	fi
	# Pulls QC count info from counts file (Order is as follows Q20_Total_[bp]	Q30_Total_[bp]	Q20_R1_[bp]	Q20_R2_[bp]	Q20_R1_[%]	Q20_R2_[%]	Q30_R1_[bp]	Q30_R2_[bp]
	# Q30_R1_[%]	Q30_R2_[%]	Total_Sequenced_[bp]	Total_Sequenced_[reads]
	read_qc_info="N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A"
	# If the counts file exists take the header line (the only one) and copy all but the first entry (which is the sample name) and store in an array
	if [[ -s "${SAMPDATADIR}/preQCcounts/${isolate_name}_counts.txt" ]]; then
		line=$(tail -n 1 "${SAMPDATADIR}/preQCcounts/${isolate_name}_counts.txt")
		IFS='	' read -r -a qcs <<< "${line}"
		read_qc_info=${qcs[@]:1}
	fi

	source_call=$(head -n1 "${SAMPDATADIR}/${isolate_name}.tax")
	tax_source="UNK"
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line:0:1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]; then
			dec_species=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "G" ]; then
			dec_genus=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "(" ]; then
			tax_source=$(echo "${line}" | cut -d'(' -f2 | cut -d')' -f1)
		fi
	done < "${SAMPDATADIR}/${isolate_name}.tax"

	# Pulls contig info from toms qc analysis file
	contig_info="0(0)\\t0\tNot_in_DB"
	if [[ -s "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv" ]]; then
		N50=$(sed -n '18p' "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv"  | sed -r 's/[\t]+/ /g'| cut -d' ' -f2)
		num_contigs=$(sed -n '14p' "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv"| sed -r 's/[\t]+/ /g' | cut -d' ' -f3 )
		assembly_length=$(sed -n '16p' "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
		#Check Assembly ratio against expected size to see if it is missing a large portion or if there is contamination/double genome
		dec_genus_initial="${dec_genus:0:1}"
		if [[ "${dec_genus_initial}" = "[" ]] || [[ "${dec_genus_initial}" = "(" ]]; then
			dec_genus_initial="${dec_genus:1:1}"
		fi
		assembly_ID="${dec_genus_initial}.${dec_species}"
		# echo "About to check mmb_bugs[${assembly_ID}]"
		# if [[ ! -z "${mmb_bugs[${assembly_ID}]}" ]]; then
		# 	assembly_ratio=$(awk -v p="${assembly_length}" -v q="${mmb_bugs[${assembly_ID}]}" -v r="${tax_source}" -v s="${assembly_ID}" 'BEGIN{printf("%.2f(%s-%s)",p/q, r, s)}')
		# else
		# 	assembly_ratio="Not_in_DB (${tax_source}-${assembly_ID})"
		# fi
		echo -e "Checking if directories exist:\nParent:${SAMPDATADIR}\nANI:${SAMPDATADIR}/ANI\nAssembly:${SAMPDATADIR}/Assembly"
		# Checks for proper argumentation
		if [[ ! -d "${SAMPDATADIR}/ANI" ]] || [[ ! -d "${SAMPDATADIR}/Assembly" ]]; then
			echo "No sample (or ANI or Assembly) folder exist, exiting"
		# Set to something to note which folder is missing?
		else
			echo "Checking if Assembly_stats exists:${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv"
			while IFS='' read -r line; do
				IFS=$'\t' read -a arr_line <<< "$line"
				#echo "${arr_line[0]}"
				#echo  "${genus} ${species} vs ${arr_line[0]}"
				if [[ "${dec_genus} ${dec_species}" = "${arr_line[0]}" ]]; then
					taxid="${arr_line[19]}"
					if [ "${taxid}" = -2 ]; then
						taxid="No mode available when determining tax id"
					elif [ "${taxid}" = -1 ]; then
						taxid="No tax id given or empty when making lookup"
					fi
					expected_length=$(echo "scale=0; 1000000 * ${arr_line[4]} / 1 " | bc | cut -d'.' -f1)
					#echo "${arr_line[5]}"
					stdev=$(echo "scale=4; 1000000 * ${arr_line[5]} /1 " | bc | cut -d"." -f1)
					if [[ "${stdev}" = "0" ]]; then
						stdev="Single_Reference"
						stdevs=0
					else
						if [[ "${assembly_length}" -gt "${expected_length}" ]]; then
							bigger="${assembly_length}"
							smaller="${expected_length}"
						else
							smaller="${assembly_length}"
							bigger="${expected_length}"
						fi
						stdevs=$(echo "scale=4 ; ( ${bigger} - ${smaller} ) / ${stdev}" | bc)
					fi
					break
				elif [[ "${dec_genus:0:1}" < "${arr_line[0]:0:1}" ]]; then
					# Cut losses if Genus initial is larger than reference list initial
					break
				fi
			done < "${NCBI_ratio}"
			#echo "looked in ${NCBI_ratio}"

			if [[ ! ${expected_length} ]]; then
				echo "No expected length was found to compare to"
				echo -e "Tax: ${dec_genus} ${dec_species}\nNCBI_TAXID: ${taxid}\nSpecies_StDev: ${stdev}\nIsolate_St.Devs: ${stdevs}\nActual_length: ${assembly_length}\nExpected_length: ${expected_length}\nRatio: -1" >  "${SAMPDATADIR}/${isolate_name}_Assembly_ratio_${NCBI_ratio_date}.txt"
			elif [[ ! ${assembly_length} ]]; then
				echo "No assembly length was found to compare with"
				echo -e "Tax: ${dec_genus} ${dec_species}\nNCBI_TAXID: ${taxid}\nSpecies_StDev: ${stdev}\nIsolate_St.Devs: ${stdevs}\nActual_length: ${assembly_length}\nExpected_length: ${expected_length}\nRatio: -2" >  "${SAMPDATADIR}/${isolate_name}_Assembly_ratio_${NCBI_ratio_date}.txt"
			fi

			ratio=$(echo "scale=6; ${assembly_length} / ${expected_length}" | bc | awk '{printf "%.4f", $0}')

			echo -e "Actual - ${assembly_length}\nExpected - ${expected_length}\nRatio - ${ratio}\nSpecies_St.Devs - ${stdev}\nIsolate_St.Dev:${stdevs}"
			echo -e "Tax: ${dec_genus} ${dec_species}\nNCBI_TAXID: ${taxid}\nSpecies_St.Dev: ${stdev}\nIsolate_St.Devs: ${stdevs}\nActual_length: ${assembly_length}\nExpected_length: ${expected_length}\nRatio: ${ratio}" >  "${SAMPDATADIR}/${isolate_name}_Assembly_ratio_${NCBI_ratio_date}.txt"
		fi
	else
		echo "No Assembly_Stats exists, so cannot continue calculate assembly ratio"
	fi
	# Set to something to note assembly stats is missing
	contig_info=$(echo -e "${num_contigs}\\t${assembly_length}\\t${assembly_ratio}")

	# Pulls busco info from summary file
	busco_info="No BUSCO performed"
	if [[ -s "${SAMPDATADIR}/BUSCO/short_summary_${isolate_name}.txt" ]]; then
		while IFS= read -r line; do
			if [[ ${line} == *"Complete BUSCOs (C)"* ]]
			then
				#echo "C-"${line}
				found_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
			elif [[ ${line} == *"Total BUSCO groups searched"* ]];
			then
				#echo "T-"${line}
				total_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
			elif [[ "${line}" == *"The lineage dataset is:"* ]];
			then
				#echo "L-"${line}
				db=$(echo "${line}" | awk -F ' ' '{print $6}')
			fi
		done < ${SAMPDATADIR}/BUSCO/short_summary_${isolate_name}.txt
		busco_info="${found_buscos}/${total_buscos}(${db})"
	fi
	# Pulls ANI info from best_ANI_hits file
	ani_info="No ANI performed"
	# Count the number of matching format files for the current sample
	file_count=$(find "${SAMPDATADIR}/ANI/" -name *"${isolate_name}"*"_vs_"*".txt" | wc -l)

	# If 1 and only 1 file exists pull the first line as the best hit information
	if [[ -s "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_REFSEQ_${REFSEQ_date}.txt" ]]; then
		ani_info=$(head -n 1 "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_REFSEQ_${REFSEQ_date}).txt")
	elif [[ -s "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_${dec_genus}).txt" ]]; then
		ani_info=$(head -n 1 "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_${dec_genus}).txt")
	# Report that more than one file exists
	else
		for file in "${SAMPDATADIR}/ANI/"*
		do
			if [[ "${file}" == *"best_ANI_hits_ordered(${isolate_name}_vs_"* ]]; then
				filename=${file}
				ani_info=$(head -n 1 "${file}")
				break
			fi
		done
	fi
	# Extract q30 reads from qcCounts to calculate average coverage as q30_reads/assembly_length
	q30_reads=$(echo "${read_qc_info}" | awk -F ' ' '{print $2}')
	# Change later to AWK as this wont work on ASPEN, but consolidate won't likely be run on cluster
	if [[ ${assembly_length} -gt 0 ]]; then
		avg_coverage=$(bc <<<"scale=2 ; ${q30_reads} / ${assembly_length}")
	else
		avg_coverage="N/A"
	fi
	# Replace all spaces in qc_info as tabs
	read_qc_info=$(echo "${read_qc_info}" | tr '[:blank:]' \\t)

	# Get the date to show when the log was made
	NOW=$(date +"%m/%d/%Y")

	# Add all pertinent info to the output file in the correct formatting to add to MMB_Seq log
	echo -e "${isolate_name}\\t${NOW}\\t${g_s_reads}\\t${g_s_assembled}\\t${g_s_16s}\\t${read_qc_info}\\t${avg_coverage}\\t${contig_info}\\t${busco_info}\\t${ani_info}\\r" >> "${PROJDATADIR}/Seqlog_output.txt"
done < ${list_path}

# Task: Create run summary
write_Progress
run_task_id=10
runsumdate=$(date "+%m_%d_%Y_at_%Hh_%Mm")
${src}/run_sum.sh ${PROJDATADIR} ${src} ${local_DBs}

# Task: Create matrix summary
write_Progress
run_task_id=11
${src}/run_MATRIX.sh -r ${PROJECT} -c ${config_file}

# Task: Copy config file to run folder to show configuration used in the run
write_Progress
run_task_id=12
echo "Moving config file(${config_file}) to log directory ($log_dir)"
cp "${config_file}" "${log_dir}/config_${PROJECT}.sh"



write_Progress
run_task_id=13
end_date=$(date "+%m_%d_%Y_at_%Hh_%Mm")
echo "Run ended at ${end_date}" #>> "${log_file}"

conda deactivate
