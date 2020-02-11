#!/bin/bash -l

#$ -o quaisar_X.out
#$ -e quaisar_X.err
#$ -N quasX
#$ -cwd
#$ -q short.q

#
# Description: The full QuAISAR-H pipeline start to end serially
#
# Usage: ./quaisar_containerized.sh -c location_of_config_file -i location_of_reads 1|2|3|4 -o path_to_parent_output_folder_location name_of_output_folder [-a]"
#		filename postfix numbers are as follows 1:_L001_SX_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz"
#
# Output location: default_config.sh_output_location
#
# Modules required: Python3/3.5.2
#
# v1.0 (2/6/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checking for proper number of arguments from command line
if [[ $# -lt 1  || $# -gt 8 ]]; then
	echo "If reads are in default location set in config file then"
  echo "Usage: ./quaisar_containerized.sh -c location_of_config_file -i location_of_reads 1|2|3|4 -o path_to_parent_output_folder_location name_of_output_folder [-a]"
	echo "filename postfix numbers are as follows 1:_L001_SX_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz"
  echo "You have used $# args"
  exit 3
fi

# Checks the arguments (more to come)
nopts=$#
do_download="false"
global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
requestor=$(whoami)
PROJECT="${requestor}_${global_time}"
assemblies="false"

for ((i=1 ; i <= nopts ; i++)); do
	#echo "${1} ${2}"
	case "${1}" in
		#Help/Usage section
		-h | --help)
			echo -e "\\n\\n\\n"
			echo "If reads are in default location set in config file then"
		  echo "Usage: ./parallel_quaisar.sh -p project_name"
			echo "else if you are running it on reads not in the default location or format"
			echo "Usage: ./quaisar_containerized.sh -c location_of_config_file -i location_of_reads 1|2|3|4 -o path_to_parent_output_folder_location name_of_output_folder [-a]"
			echo "filename postfix numbers are as follows 1:_L001_SX_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz"
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
			postfix_index="$3"
			#is_full_run="false"
			#echo "$INDATADIR $2"
			shift 3
			;;
		#Gets output directory name of folder that all output files will be stored
		-o | --out-dir)
			BASEDIR="$2"
			PROJECT="$3"
			shift 3

			echo "output_dir=${BASEDIR}" >> "${src}/config.sh"
			. ${src}/config.sh
			echo "${output_dir}"
			list_path="${BASEDIR}/${PROJECT}/${PROJECT}_list.txt"
			if [[ ! -d ${BASEDIR} ]]; then
				mkdir -p ${BASEDIR}
			fi
			;;
		-a | --assemblies)
				assemblies="true"
				shift
				;;
		#Captures any other characters in the args
		\?)
			echo "ERROR: ${BOLD}$2${NORM} is not a valid argument" >&2
			usage
			exit 1
			;;
	esac
done

# Short print out summary of run settings
echo -e "Source folder: ${INDATADIR}\\nOutput folder: ${BASEDIR}\\nList based analysis:  ${list_path}"

# Sets folder to where files will be downloaded to
PROJDATADIR="${output_dir}/${PROJECT}"
if [ ! -d "${PROJDATADIR}" ]; then
	echo "Creating $PROJDATADIR"
	mkdir -p "${PROJDATADIR}"
fi
if [ -f "${PROJDATADIR}/${PROJECT}_list.txt" ]; then
	mv "${PROJDATADIR}/${PROJECT}_list.txt" "${PROJDATADIR}/${PROJECT}_list_original.txt"
fi


# Copies reads from source location to working directory and creates a list of IDs
if [[ "${assemblies}" == "true" ]]; then
	# Goes through given Assemblies folder
	echo "${INDATADIR}"
	for file in ${INDATADIR}/*
	do
		# Check if file is a recognized assembly format externsion
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
		if [[ "${file}" = *.gz ]] || [[ "${file}" = *.fastq ]]; then
			echo "isolate_name: ${file}"
			# Gets full file name from path
			full_sample_name=${file##*/}
			if [[ "${full_sample_name}" = *"_R1_"* ]]; then
				full_sample_name_pair=${full_sample_name/_R1_/_R2_}
			elif [[ "${full_sample_name}" = *"_R2_"* ]]; then
				full_sample_name_pair="${full_sample_name/_R2_/_R1_}"
			elif [[ "${full_sample_name}" = *"_1.fast"* ]]; then
				full_sample_name_pair=${full_sample_name/_1.fast/_2.fast}
			elif [[ "${full_sample_name}" = *"_2.fast"* ]]; then
				full_sample_name_pair="${full_sample_name/_2.fast/_1.fast}"
			fi
			# gets path from file
			source_path=$(dirname "${file}")
			# Extracts isolate_name keeping only isolate ID, if it matches standard miseq naming
			echo "${postfix_index}:${full_sample_name}"
			if [[ "${postfix_index}" -eq 1 ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f5- | rev)
				postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2,3 | rev)
			elif [[ "${postfix_index}" -eq 4 ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f4- | rev)
				postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2,3 | rev)
			elif [[ "${postfix_index}" -eq 3 ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f3- | rev)
				postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2 | rev)
			elif [[ "${postfix_index}" -eq 2 ]]; then
				short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f2- | rev)
				postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1 | rev)
			else
				echo "Magic - should have never gotten here as this number does not match any of the input numbers... 1:_L001_SX_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz , exiting"
				exit
			fi

			#long_name=$(echo "${full_sample_name}" | cut -d'_' -f1,2,3)
			echo "Short: ${short_name}"
			#echo "Does ${full_sample_name} match *${match}"

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
					echo "Creating $PROJDATADIR/${short_name}/FASTQs"
					mkdir -p "${PROJDATADIR}/${short_name}/FASTQs"
				fi
				# Announces name of file being unzipped and then unzips it to the FASTQs folder for the matching sample name. Files are shortened to just name_R1_001.fastq or name_R2_001.fastq
				echo "Retrieving ${source_path}/${full_sample_name} and ${full_sample_name_pair}"
				#if [[ "${match}" -eq 1 ]] || [[ "${match}" -eq 4 ]]; then
					if [[ "${postfix}" = *"R1_001.fast"* ]] || [[ "${postfix}" = *"R1.fast"* ]] || [[ "${postfix}" = *"1.fast"* ]]; then
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
					fi
					if [[ "${postfix}" = *"R2_001.fast"* ]] || [[ "${postfix}" = *"R2.fast"* ]] || [[ "${postfix}" = *"2.fast"* ]]; then
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
						fi
					fi
					if grep -Fxq "${PROJECT}/${short_name}" "${PROJDATADIR}/${PROJECT}_list.txt"
					then
						echo -e "${PROJECT}/${short_name} already on list "${PROJDATADIR}/${PROJECT}_list.txt", not adding again"
					else
						echo -e "${PROJECT}/${short_name}" >> "${PROJDATADIR}/${PROJECT}_list.txt"
					fi
			fi
		else
			echo "${file} is not a FASTQ(.gz) read file, not acting on it"
		fi
	done
fi

# Invert list so that the important isolates (for us at least) get run first
if [[ -f "${PROJDATADIR}/${PROJECT}_list.txt" ]]; then
	sort -k2,2 -t'/' -r "${PROJDATADIR}/${PROJECT}_list.txt" -o "${PROJDATADIR}/${PROJECT}_list.txt"
fi

# Loops through list file to create an array of all isolates to run through pipeline
declare -a isolate_list=()
while IFS= read -r file || [ -n "$file" ]; do
	echo "Found: ${file}"
	file=$(echo "${file}" | tr -d '[:space:]')
	isolate_list+=("${file}")
done < "${list_path}"

# Displays number and names of files found to analyze
if [[ ${#isolate_list[@]} -gt 1 ]]; then
	echo "Will analyze these ${#isolate_list[@]} files: ${isolate_list[*]}"
elif [[ ${#isolate_list[@]} -eq 1 ]]; then
	echo "Will analyze this file: ${isolate_list[0]}"
else
	echo "No files found in ${list_path}"
fi

run_start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

#Each file in the list is checked individually for successful completion and added then added to the log for the run
log_dir="${PROJDATADIR}"

#Get the time the run started to use as the identifier
outarray=()
echo "Run started at ${run_start_time}; Log saved to ${log_dir}/${PROJECT}_on_${run_start_time}.log"
echo "Run started at ${run_start_time}" > "${log_dir}/${PROJECT}_on_${run_start_time}.log"
outarray+=("${PROJECT} started at ${run_start_time} and saved to ${log_dir}/${PROJECT}_on_${run_start_time}.log")

for isolate in "${isolate_list[@]}"; do

	#Time tracker to gauge time used by each step
	totaltime=0
	start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

	# Set arguments to isolate_name PROJECT (miseq run id) and SAMPDATADIR(${output_dir}/PROJECT/isolate_name)
	isolate_name=$(echo "${isolate}" | awk -F/ '{print $2}' | tr -d '[:space:]')
	SAMPDATADIR="${PROJDATADIR}/${isolate_name}"

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
		#Checks if FASTQ folder exists for current sample
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
		# Checks for and creates the specified output folder for the QC counts
		if [ ! -d "${SAMPDATADIR}/preQCcounts" ]; then
			echo "Creating ${SAMPDATADIR}/preQCcounts"
			mkdir -p "${SAMPDATADIR}/preQCcounts"
		fi
		# Run qc count check on raw reads
		python3 "${src}/Fastq_Quality_Printer.py" -1 "${SAMPDATADIR}/FASTQs/${isolate_name}_R1_001.fastq" -2 "${SAMPDATADIR}/FASTQs/${isolate_name}_R2_001.fastq" > "${SAMPDATADIR}/preQCcounts/${isolate_name}_counts.txt"
		# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeQCcount=$((end - start))
		echo "QC count - ${timeQCcount} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeQCcount))

		###  Trimming and Quality Control  ###
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
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES ${src}/singularity_images/bbtools.simg bbduk.sh -${bbduk_mem} threads=${procs} in=/SAMPDIR/FASTQs/${isolate_name}_R1_001.fastq in2=/SAMPDIR/FASTQs/${isolate_name}_R2_001.fastq out=/SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R1.fsq out2=/SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R2.fsq ref=/DATABASES/phiX.fasta k=${bbduk_k} hdist=${bbduk_hdist}
		# Get end time of bbduk and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeAdapt=$((end - start))
		echo "Removing Adapters - ${timeAdapt} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeAdapt))


		### Quality and Adapter Trimming using trimmomatic ###
		echo "----- Running Trimmomatic on reads -----"
		# Get start time of trimmomatic
		start=$SECONDS
		# Creates folder for trimmomatic output if it does not exist
		if [ ! -d "${SAMPDATADIR}/trimmed" ]; then
			mkdir -p "${SAMPDATADIR}/trimmed"
		fi
		### Run trimmomatic
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/trimmomatic:0.36--5 trimmomatic ${trim_endtype} -${trim_phred} -threads ${procs} /SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R1.fsq /SAMPDIR/removedAdapters/${isolate_name}-noPhiX-R2.fsq /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq /SAMPDIR/trimmed/${isolate_name}_R1_001.unpaired.fq /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq /SAMPDIR/trimmed/${isolate_name}_R2_001.unpaired.fq ILLUMINACLIP:/DATABASES/adapters.fasta:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome} SLIDINGWINDOW:${trim_window_size}:${trim_window_qual} LEADING:${trim_leading} TRAILING:${trim_trailing} MINLEN:${trim_min_length}
		# Get end time of trimmomatic and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeTrim=$((end - start))
		echo "Trimming - ${timeTrim} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeTrim))
		gzip -k "${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq"
		gzip -k "${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq"

		# Check differences after QC and trimming (also for gottcha proper read count for assessing unclassified reads)
		# Get start time for qc check on trimmed reads
		start=$SECONDS
		### Count the number of Q20, Q30, bases and reads within the trimmed pair of FASTQ files
		echo "----- Counting read quality of trimmed files-----"
		# Run qc count check on filtered reads
		python3 "${src}/Fastq_Quality_Printer.py" -1 "${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq" -2 "${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq" > "${SAMPDATADIR}/preQCcounts/${isolate_name}_trimmed_counts.txt"
		# Merge both unpaired fq files into one for GOTTCHA
		cat "${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.unpaired.fq" "${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.unpaired.fq" > "${SAMPDATADIR}/trimmed/${isolate_name}.single.fq"
		cat "${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq" "${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq" > "${SAMPDATADIR}/trimmed/${isolate_name}.paired.fq"

		# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeQCcount_trimmed=$((end - start))
		echo "QC count trimmed - ${timeQCcount_trimmed} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeQCcount))

		######  Run Kraken on cleaned reads  ######
		echo "----- Running Kraken on cleaned reads -----"
		# Get start time of kraken on reads
		start=$SECONDS
		# Create directory for kraken dataset
		if [[ ! -d "${SAMPDATADIR}/kraken/preAssembly" ]]; then
			mkdir -p "${SAMPDATADIR}/kraken/preAssembly"
		fi
		# Run kraken
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken --paired --db /DATABASES/minikrakenDB/  --preload --fastq-input --threads 4 --output /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.kraken --classified-out /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.classified /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken-mpa-report --db /DATABASES/minikrakenDB/ /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.kraken > "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.mpa"
		python3 "${src}/Metaphlan2krona.py" -p "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.mpa" -k "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.krona"
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR docker://quay.io/biocontainers/krona:2.7--0 ktImportText /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.krona -o /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.html
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken-report --db /DATABASES/minikrakenDB/ /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.kraken > "${SAMPDATADIR}/kraken/preAssembly/${isolate_name}_paired.list"
		"${src}/best_hit_from_kraken.sh" "${isolate_name}" "pre" "paired" "${PROJECT}" "kraken"
		# Get end time of kraken on reads and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeKrak=$((end - start))
		echo "Kraken - ${timeKrak} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeKrak))

		##### Run gottcha(v1) on cleaned reads #####
		echo "----- Running gottcha on cleaned reads -----"
		# Get start time of gottcha
		start=$SECONDS
		# run gottcha
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B "${local_DBs}":/DBs ${src}/singularity_images/gottcha.simg gottcha.pl --mode all --outdir /SAMPDIR/gottcha/gottcha_S --input /SAMPDIR/trimmed/${isolate_name}.paired.fq --database /DBs/gottcha/GOTTCHA_BACTERIA_c4937_k24_u30.species
		singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR docker://quay.io/biocontainers/krona:2.7--0 ktImportText /SAMPDIR/gottcha/gottcha_S/${isolate_name}_temp/${isolate_name}.lineage.tsv -o /SAMPDIR/gottcha/${isolate_name}_species.krona.html
		#Create a best hit from gottcha1 file
		"${src}/best_hit_from_gottcha1.sh" "${isolate_name}" "${PROJECT}"
		# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
		end=$SECONDS
		timeGott=$((end - start))
		echo "Gottcha - ${timeGott} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeGott))

		# Check reads using SRST2
		echo "----- Running SRST2 -----"
		start=$SECONDS

		if [[ ! -d "${SAMPDATADIR}/srst2" ]]; then
				mkdir "${SAMPDATADIR}/srst2"
		fi

		cp ${SAMPDATADIR}/trimmed/${isolate_name}_R1_001.paired.fq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz
		cp ${SAMPDATADIR}/trimmed/${isolate_name}_R2_001.paired.fq.gz ${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz

		singularity exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DBs ${src}/singularity_images/srst2.simg srst2 --input_pe /SAMPDIR/srst2/testreads_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/testreads_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/testreads_ResGANNCBI_20191227 --gene_db /DBs/star/ResGANNCBI_20191227_srst2.fasta

		# Cleans up leftover files
		rm "${SAMPDATADIR}/srst2/"*".bam"
		rm "${SAMPDATADIR}/srst2/"*".pileup"
		rm "${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz"
		rm "${SAMPDATADIR}/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz"

		# Removes the extra ResGANNCBI__ from all files created
		find ${SAMPDATADIR}/srst2 -type f -name "*ResGANNCBI__*" | while read FILE ; do
		  dirname=$(dirname $FILE)
			filename=$(basename $FILE)
			filename="${filename/_ResGANNCBI__/__}"
			#echo "Found-${FILE}"
			#echo "${filename}"
		    mv "${FILE}" "${dirname}/${filename}"
		done
		end=$SECONDS
		timesrst2=$((end - start))
		echo "SRST2 - ${timesrst2} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timesrst2))

		######  Assembling Using SPAdes  ######
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
				singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/spades:3.13.0--0 spades.py --careful --memory 16 --only-assembler --pe1-1 /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq --pe1-2 /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq --pe1-s /SAMPDIR/trimmed/${isolate_name}.single.fq -o /SAMPDIR/Assembly --phred-offset "${phred}" -t "${procs}"
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
	fi

	### Removing Short Contigs  ###
	echo "----- Removing Short Contigs -----"
	python3 "${src}/removeShortContigs.py" -i "${SAMPDATADIR}/Assembly/scaffolds.fasta" -t 500 -s "normal_SPAdes"
	mv "${SAMPDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta"

	# Checks to see that the trimming and renaming worked properly, returns if unsuccessful
	if [ ! -s "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" ]; then
		echo "Trimmed contigs file does not exist continuing to next sample">&2
		continue
	fi

	### ReKraken on Assembly ###
	echo "----- Running Kraken on Assembly -----"
	# Get start time of kraken on assembly
	start=$SECONDS
	# Create directory for kraken dataset on assembly
	if [[ ! -d "${SAMPDATADIR}/kraken/postAssembly" ]]; then
		mkdir -p "${SAMPDATADIR}/kraken/postAssembly"
	fi
	# Run kraken on assembly
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken --db /DATABASES/minikrakenDB/  --preload --threads ${procs} --output /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.kraken --classified-out /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.classified /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta
	python3 ${src}/Kraken_Assembly_Converter_2_Exe.py -i "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.kraken"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken-mpa-report --db /DATABASES/minikrakenDB/ /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled_BP.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted.mpa"
	python3 "${src}/Metaphlan2krona.py" -p "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted.mpa" -k "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_weighted.krona"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken-report --db /DATABASES/minikrakenDB/ /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled_BP.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled_BP.list"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/krona:2.7--0 ktImportText /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled_weighted.krona -o /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled_weighted_BP_krona.html
	"${src}/best_hit_from_kraken.sh" "${isolate_name}" "post" "assembled_BP" "${PROJECT}" "kraken"singularity -s exec -B ${SAMPDATADIR}/kraken/postAssembly:/INPUT -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken-report --db /DATABASES/minikrakenDB/ /INPUT/${isolate_name}_assembled.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.list"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken-mpa-report --db /DATABASES/minikrakenDB/ /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.mpa"
	python3 "${src}/Metaphlan2krona.py" -p "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.mpa" -k "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.krona"
	singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR docker://quay.io/biocontainers/krona:2.7--0 ktImportText /SAMPDIR/kraken/preAssembly/${isolate_name}_paired.krona -o /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.html
	singularity -s exec -B "${SAMPDATADIR}":/SAMPDIR -B ${local_DBs}:/DATABASES docker://quay.io/biocontainers/kraken:1.0--pl5.22.0_0 kraken-report --db /DATABASES/minikrakenDB/ /SAMPDIR/kraken/postAssembly/${isolate_name}_assembled.kraken > "${SAMPDATADIR}/kraken/postAssembly/${isolate_name}_assembled.list"
	"${src}/best_hit_from_kraken.sh" "${isolate_name}" "post" "assembled" "${PROJECT}" "kraken"
	# Get end time of kraken on assembly and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeKrakAss=$((end - start))
	echo "Kraken on Assembly - ${timeKrakAss} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeKrakAss))

	# Get ID fom 16s
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
	#barrnap --kingdom bac --threads ${procs} "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" > ${SAMPDATADIR}/16s/${isolate_name}_scaffolds_trimmed.fasta_rRNA_seqs.fasta

	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/barrnap:0.8--0 barrnap --kingdom bac /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta > "${SAMPDATADIR}/16s/${isolate_name}_scaffolds_trimmed.fasta_rRNA_seqs.fasta"

	# Checks for successful output from barrnap, *rRNA_seqs.fasta
	if [[ ! -s ${SAMPDATADIR}/16s/${isolate_name}_scaffolds_trimmed.fasta_rRNA_seqs.fasta ]]; then
		echo "rNA_seqs.fasta does NOT exist"
	fi

	# Checks barrnap output and finds all 16s hits and creates a multi-fasta file to list all possible matches
	lines=0
	found_16s="false"
	while IFS='' read -r line || [ -n "$line" ]; do
		if [ ${lines} -gt 0 ]; then
			contig=$(echo ${line} | cut -d' ' -f1)
			cstart=$(echo ${line} | cut -d' ' -f4)
			cstop=$(echo ${line} | cut -d' ' -f5)
			ribosome=$(echo ${line} | cut -d' ' -f9 | cut -d'=' -f3)
			if [ "${ribosome}" = "16S" ]; then
				# Replace with subsequence once it can handle multi-fastas
				#make_fasta $1 $2 $contig $cstart $cstop
				python3 ${src}/get_subsequence.py -i "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" -s ${cstart} -e ${cstop} -t ${contig} >> ${SAMPDATADIR}/16s/${isolate_name}_16s_rna_seqs.txt
				found_16s="true"
			fi
		fi
		lines=$((lines + 1))
	done < "${SAMPDATADIR}/16s/${isolate_name}_scaffolds_trimmed.fasta_rRNA_seqs.fasta"

	# Adds No hits found to output file in the case where no 16s ribosomal sequences were found
	if [[ "${found_16s}" == "false" ]]; then
		echo -e "best_hit	${isolate_name}	No_16s_sequences_found" > "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt"
		echo -e "largest_hit	${isolate_name}	No_16s_sequences_found" >> "${SAMPDATADIR}/16s/${isolate_name}_16s_blast_id.txt"
	fi

	# Blasts the NCBI database to find the closest hit to every entry in the 16s fasta list
	###### MAX_TARGET_SEQS POSSIBLE ERROR
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/blast:2.9.0--pl526h3066fca_4 blastn -word_size 10 -task blastn -remote -db nt -max_hsps 1 -max_target_seqs 1 -query /SAMPDIR/16s/${isolate_name}_16s_rna_seqs.txt -out /SAMPDIR/16s/${isolate_name}.nt.RemoteBLASTN -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ssciname";
	# Sorts the list based on sequence match length to find the largest hit
	sort -k4 -n "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN" --reverse > "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN.sorted"

	# Gets taxon info from the best bitscore (literal top) hit from the blast list
	if [[ -s "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN" ]]; then
		me=$(whoami)
		accessions=$(head -n 1 "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN")
		hits=$(echo "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN" | wc -l)
	#	echo ${accessions}
		gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
		echo ${gb_acc}
		attempts=0
		# Will try getting info from entrez up to 5 times, as it has a higher chance of not finishing correctly on the first try
		while [[ ${attempts} -lt 5 ]]; do
			blast_id=$(singularity exec ${src}/singularity_images/entrez_taxon.simg python3 /entrez/entrez_get_taxon_from_accession.py -a "${gb_acc}" -e "${runner}")
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
		if [[ -s "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN.sorted" ]]; then
			me=$(whoami)
			accessions=$(head -n 1 "${SAMPDATADIR}/16s/${isolate_name}.nt.RemoteBLASTN.sorted")
			gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
			attempts=0
			# Will try getting info from entrez up to 5 times, as it has a higher chance of not finishing correctly on the first try
			while [[ ${attempts} -lt 5 ]]; do
				blast_id=$(singularity exec ${src}/singularity_images/entrez_taxon.simg python3 /entrez/entrez_get_taxon_from_accession.py -a "${gb_acc}" -e "${runner}")
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

	# Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
	"${src}/determine_taxID.sh" "${isolate_name}" "${PROJECT}"
	# Capture the anticipated taxonomy of the sample using kraken on assembly output
	echo "----- Extracting Taxonomy from Taxon Summary -----"
	# Checks to see if the kraken on assembly completed successfully
	if [ -s "${SAMPDATADIR}/${isolate_name}.tax" ]; then
		# Read each line of the kraken summary file and pull out each level  taxonomic unit and store for use later in busco and ANI
		while IFS= read -r line  || [ -n "$line" ]; do
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

	### Check quality of Assembly ###
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
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/QUAST5.simg python3 /quast/quast.py -o /SAMPDIR/Assembly_Stats /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta
	mv "${SAMPDATADIR}/Assembly_Stats/report.txt" "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.txt"
	mv "${SAMPDATADIR}/Assembly_Stats/report.tsv" "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv"
	# Get end time of qc quality check and calculate run time and append to time summary (and sum to total time used)
	end=$SECONDS
	timeQCcheck=$((end - start))
	echo "QC Check of Assembly - ${timeQCcheck} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeQCcheck))

	### Prokka on assembly ###
	echo "----- Running Prokka on Assembly -----"
	# Get start time for prokka
	start=$SECONDS
	# Run prokka
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/prokka:1.14.5--pl526_0 prokka --outdir /SAMPDIR/prokka /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta
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

	# Rename contigs to something helpful (Had to wait until after prokka runs due to the strict naming requirements
	mv "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed_original.fasta"
	python3 "${src}/fasta_headers.py" -i "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed_original.fasta" -o "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta"

	### Average Nucleotide Identity ###
	echo "----- Running ANI for Species confirmation -----"
	# ANI uses assembly and sample would have exited already if assembly did not complete, so no need to check
	# Get start time of ANI
	start=$SECONDS
	# run ANI
	# Temp fix for strange genera until we do vs ALL all the time.
	if [[ "${genus}" = "Peptoclostridium" ]] || [[ "${genus}" = "Clostridioides" ]]; then
		genus="Clostridium"
	elif [[ "${genus}" = "Shigella" ]]; then
		genus="Escherichia"
	fi

	if [ ! -d "${SAMPDATADIR}/ANI" ]; then  #create outdir if absent
		echo "${SAMPDATADIR}/ANI"
		mkdir -p "${SAMPDATADIR}/ANI"
	fi

	if [ ! -s "${local_DBs}/aniDB/${genus,}" ]; then
		echo "The genus does not exist in the ANI database. This will be noted and the curator of the database will be notified. However, since nothing can be done at the moment....exiting"
		# Write non-results to a file in ANI folder
		echo "No matching ANI database found for ${genus}" >> "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_${genus}).txt"
		global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
		echo "ANI: ${genus} - Found as ${isolate_name} on ${global_time}" >> "${src}/maintenance_To_Do.txt"
	fi

	# Checks to see if the local DB ANI folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
	if [ ! -d "${SAMPDATADIR}/ANI/localANIDB" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/ANI/localANIDB"
		mkdir -p "${SAMPDATADIR}/ANI/localANIDB"
	else
		rm -r "${SAMPDATADIR}/ANI/localANIDB"
		mkdir -p "${SAMPDATADIR}/ANI/localANIDB"
	fi

	# Checks to see if the local DB ANI temp folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
	if [ ! -d "${SAMPDATADIR}/ANI/temp" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/ANI/temp"
		mkdir -p "${SAMPDATADIR}/ANI/temp"
	else
		rm -r "${SAMPDATADIR}/ANI/temp"
		mkdir -p "${SAMPDATADIR}/ANI/temp"
	fi

	#Creates a local copy of the database folder
	echo "trying to copy ${local_DBs}/aniDB/${genus,}/"
	cp "${local_DBs}/aniDB/${genus,}/"*".fna.gz" "${SAMPDATADIR}/ANI/localANIDB/"
	gunzip ${SAMPDATADIR}/ANI/localANIDB/*.gz

	#Copies the samples assembly contigs to the local ANI db folder
	cp "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" "${SAMPDATADIR}/ANI/localANIDB/sample_${genus}_${species}.fasta"

	#Renames all files in the localANIDB folder by changing extension from fna to fasta (which pyani needs)
	for file in ${SAMPDATADIR}/ANI/localANIDB/*.fna;
	do
		fasta_name=$(basename "${file}" .fna)".fasta"
		mv "${file}" "${SAMPDATADIR}/ANI/localANIDB/${fasta_name}"
	done

	# Mashtree trimming to reduce run time for ANI
	owd=$(pwd)
	cd ${SAMPDATADIR}/ANI/localANIDB/

	echo "----- Running MASHTREE for inside ANI -----"
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/mashtree:1.0.4--pl526h516909a_0 mashtree --numcpus ${procs} /SAMPDIR/ANI/localANIDB/*.fasta > ${SAMPDATADIR}/ANI/"${genus}_and_${isolate_name}_mashtree.dnd"

	# Get total number of isolates compared in tree
	sample_count=$(find ${SAMPDATADIR}/ANI/localANIDB/ -type f | wc -l)
	# Must remove sample of interest
	sample_count=$(( sample_count - 1 ))
	# Check if sample count is greater than the max samples for tree size, if so then reduce tree size to max closest samples balanced around submitted isolate
	if [[ ${sample_count} -gt ${max_ani_samples} ]]; then
		sleep 2
		tree=$(head -n 1 "${SAMPDATADIR}/ANI/${genus}_and_${isolate_name}_mashtree.dnd")
		echo $tree
		tree=$(echo "${tree}" | tr -d '()')
		echo $tree
		IFS=',' read -r -a samples <<< "${tree}"
		counter=0
		half_max=$(( (max_ani_samples+1) / 2 ))
		echo "Halfsies = ${half_max}"
		for sample in ${samples[@]};
		do
			counter=$(( counter + 1 ))
			filename=$(echo ${sample} | cut -d':' -f1)
			echo "${filename}"
			filename="${filename}.fasta"
			if [[ "${filename}" == "sample_${genus}_${species}.fasta" ]]; then
				match=${counter}
				#echo "Match @ ${counter} and half=${half_max}"
				if [[ ${match} -le ${half_max} ]]; then
					#echo "LE"
					samples_trimmed=${samples[@]:0:$(( max_ani_samples + 1 ))}
				elif [[ ${match} -ge $(( sample_count - half_max)) ]]; then
					#echo "GE"
					samples_trimmed=${samples[@]:$(( max_ani_samples * -1 - 1 ))}
				else
					#echo "MID - $(( match - half_max )) to $(( counter + half_max + 1))"
					samples_trimmed=${samples[@]:$(( match - half_max )):${max_ani_samples}}
				fi
					#echo "${#samples_trimmed[@]}-${samples_trimmed[@]}"
					break
			fi
					#echo ${filename}
		done
		mkdir "${SAMPDATADIR}/ANI/localANIDB_trimmed"
		for sample in ${samples_trimmed[@]};
		do
			filename=$(echo ${sample} | cut -d':' -f1)
			filename="${filename}.fasta"
			echo "Moving ${filename}"
			cp ${SAMPDATADIR}/ANI/localANIDB/${filename} ${SAMPDATADIR}/ANI/localANIDB_trimmed/
		done
		if [[ -d "${SAMPDATADIR}/ANI/localANIDB_full" ]]; then
			rm -r "${SAMPDATADIR}/ANI/localANIDB_full"
		fi
		mv "${SAMPDATADIR}/ANI/localANIDB" "${SAMPDATADIR}/ANI/localANIDB_full"
		mv "${SAMPDATADIR}/ANI/localANIDB_trimmed" "${SAMPDATADIR}/ANI/localANIDB"
		#rm -r "${SAMPDATADIR}/ANI/localANIDB_full"

	# Continue without reducing the tree, as there are not enough samples to require reduction
	else
		echo "Sample count below limit, not trimming ANI database"
	fi

	cd ${owd}
	# Resume normal ANI analysis after mashtree reduction

	# Checks for a previous copy of the aniM folder, removes it if found
	if [ -d "${SAMPDATADIR}/ANI/aniM" ]; then
		echo "Removing old ANIm results in ${SAMPDATADIR}/ANI/aniM"
		rm -r "${SAMPDATADIR}/ANI/aniM"
	fi

	#Calls pyani on local db folder
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/pyani:0.2.7--py35h24bf2e0_1 average_nucleotide_identity.py -i /SAMPDIR/ANI/localANIDB -o /SAMPDIR/ANI/aniM --write_excel

	#Extracts the query sample info line for percentage identity from the percent identity file
	while IFS='' read -r line; do
	#	echo "!-${line}"
		if [[ ${line:0:7} = "sample_" ]]; then
			sampleline=${line}
	#		echo "found it-"$sampleline
			break
		fi
	done < "${SAMPDATADIR}/ANI/aniM/ANIm_percentage_identity.tab"

	#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
	if [[ -s "${SAMPDATADIR}/ANI/aniM/ANIm_percentage_identity.tab" ]]; then
		firstline=$(head -n 1 "${SAMPDATADIR}/ANI/aniM/ANIm_percentage_identity.tab")
	else
		echo "No ${SAMPDATADIR}/ANI/aniM/ANIm_percentage_identity.tab file, exiting"
	fi

	#Arrays to read sample names and the %ids for the query sample against those other samples
	IFS="	" read -r -a samples <<< "${firstline}"
	IFS="	" read -r -a percents <<< "${sampleline}"

	#How many samples were compared
	n=${#samples[@]}

	#Extracts all %id against the query sample (excluding itself) and writes them to file
	for (( i=0; i<n; i++ ));
	do
	#	echo ${i}-${samples[i]}
		if [[ ${samples[i]:0:7} = "sample_" ]];
		then
	#		echo "Skipping ${i}"
			continue
		fi
		definition=$(head -1 "${SAMPDATADIR}/ANI/localANIDB/${samples[i]}.fasta")
		# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
		echo "${percents[i+1]} ${samples[i]} ${definition}" >> "${SAMPDATADIR}/ANI/best_hits.txt"
	done

	#Sorts the list in the file based on %id (best to worst)
	sort -nr -t' ' -k1 -o "${SAMPDATADIR}/ANI/best_hits_ordered.txt" "${SAMPDATADIR}/ANI/best_hits.txt"
	#Extracts the first line of the file (best hit)
	best=$(head -n 1 "${SAMPDATADIR}/ANI/best_hits_ordered.txt")
	#Creates an array from the best hit
	IFS=' ' read -r -a def_array <<< "${best}"
	#Captures the assembly file name that the best hit came from
	best_file=${def_array[1]}
	#Formats the %id to standard percentage (xx.xx%)
	best_percent=$(awk -v per="${def_array[0]}" 'BEGIN{printf "%.2f", per * 100}')
	#echo "${best_file}"
	#Extracts the accession number from the definition line
	accession=$(echo "${def_array[2]}" | cut -d' ' -f1  | cut -d'>' -f2)
	#Looks up the NCBI genus species from the accession number
	if [[ "${accession}" == "No_Accession_Number" ]]; then
		best_organism_guess="${def_array[3]} ${def_array[4]}"
	else
		attempts=0
		while [[ ${attempts} -lt 25 ]]; do
			best_organism_guess=$(singularity exec ${src}/singularity_images/entrez_taxon.simg python3 /entrez/entrez_get_taxon_from_accession.py -a "${accession}" -e "${runner}")
			if [[ ! -z ${best_organism_guess} ]]; then
				break
			else
				attempts=$(( attempts + 1 ))
			fi
		done
	fi

	#Creates a line at the top of the file to show the best match in an easily readable format that matches the style on the MMB_Seq log
	echo -e "${best_percent}%-${best_organism_guess}(${best_file}.fna)\\n$(cat "${SAMPDATADIR}/ANI/best_hits_ordered.txt")" > "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_${genus_in}).txt"

	# Get end time of ANI and calculate run time and append to time summary (and sum to total time used
	end=$SECONDS
	timeANI=$((end - start))
	echo "autoANI - ${timeANI} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeANI))

	# Get taxonomy from currently available files
	"${src}/determine_taxID.sh" "${isolate_name}" "${PROJECT}"

	### BUSCO on prokka output ###
	echo "----- Running BUSCO on Assembly -----"
	# Check to see if prokka finished successfully
	if [ -s "${SAMPDATADIR}/prokka/${isolate_name}_PROKKA.gbf" ] || [ -s "${SAMPDATADIR}/prokka/${isolate_name}_PROKKA.gff" ]; then
		# Get start time of busco
		start=$SECONDS
		# Set default busco database as bacteria in event that we dont have a database match for sample lineage
		buscoDB="bacteria_odb10"
		# Iterate through taxon levels (species to domain) and test if a match occurs to entry in database. If so, compare against it
		busco_found=0
		for tax in $species $genus $family $order $class $phylum $kingdom $domain
		do
			if [ -d "${local_DBs}/BUSCO/${tax,}_odb10" ]
			then
				buscoDB="${tax,}_odb10"
				busco_found=1
				break
			fi
		done
		# Report an unknown sample to the maintenance file to look into
		if [[ "${busco_found}" -eq 0 ]]; then
			global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
			echo "BUSCO: ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species} - Found as ${PROJECT}/${isolate_name} on ${global_time}" >> "${src}/maintenance_To_Do.txt"
		fi
		# Show which database entry will be used for comparison
		echo "buscoDB:${buscoDB}"
		# Run busco
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}/BUSCO:/DATABASES docker://quay.io/biocontainers/busco:3.0.2--py35_4 run_BUSCO.py -i /SAMPDIR/prokka/${isolate_name}_PROKKA.faa -o /SAMPDIR/BUSCO -l /DATABASES/${buscoDB} -m prot

		# Get end time of busco and calculate run time and append to time summary (and sum to total time used
		end=$SECONDS
		timeBUSCO=$((end - start))
		echo "BUSCO - ${timeBUSCO} seconds" >> "${time_summary}"
		totaltime=$((totaltime + timeBUSCO))
	# Prokka did not complete successfully and busco cant run (since it relies on prokka output)
	else
		echo "Prokka output not found, not able to process BUSCO"
	fi

	### c-SSTAR for finding AR Genes ###
	echo "----- Running c-SSTAR for AR Gene identification -----"
	# c-SSTAR uses assembly and sample would have exited already if assembly did not complete, so no need to check
	# Get start time of ccstar
	start=$SECONDS

	# Run csstar in default mode from config.sh
	if [ ! -d "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/c-sstar${ResGANNCBI_srst2_filename}_${csstar_gapping}"
		mkdir -p "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}"
	fi
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES ${src}/singularity_images/cSSTAR.simg python3 /cSSTAR/c-SSTAR_${csstar_gapping}.py -g /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -s "${csim}" -d /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta > "${SAMPDATADIR}/c-sstar/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${csim}.sstar"

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
	end=$SECONDS
	timestar=$((end - start))
	echo "c-SSTAR - ${timestar} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timestar))

	start=$SECONDS
	# Run GAMA on assembly
	echo "----- Running GAMA -----"
	if [ ! -d "${SAMPDATADIR}/GAMA" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/GAMA"
		mkdir -p "${SAMPDATADIR}/GAMA"
	fi
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR -B ${local_DBs}:/DATABASES ${src}/singularity_images/GAMA_quaisar.simg python3 /GAMA/GAMA_quaisar.py -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -d /DATABASES/star/${ResGANNCBI_srst2_filename}_srst2.fasta -o /SAMPDIR/GAMA/${isolate_name}.${ResGANNCBI_srst2_filename}.GAMA

	end=$SECONDS
	timeGAMA=$((end - start))
	echo "GAMA - ${timeGAMA} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeGAMA))

	# Get MLST profile
	echo "----- Running MLST -----"
	if [ ! -d "${SAMPDATADIR}/MLST" ]; then  #create outdir if absent
		echo "Creating ${SAMPDATADIR}/MLST"
		mkdir -p "${SAMPDATADIR}/MLST"
	fi
	start=$SECONDS
	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/mlst:2.16--0 mlst /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta > "${SAMPDATADIR}/MLST/${isolate_name}.mlst"
	python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}.mlst" -t standard
	type=$(head -n1 ${SAMPDATADIR}/MLST/${isolate_name}.mlst | cut -d' ' -f3)
	if [[ "${genus}_${species}" = "Acinetobacter_baumannii" ]]; then
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/mlst:2.16--0 mlst --scheme "abaumannii" /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta > "${SAMPDATADIR}/MLST/${isolate_name}_abaumannii.mlst"
		python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_abaumannii.mlst" -t standard
		mv "${SAMPDATADIR}/MLST/${isolate_name}_abaumannii.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Oxford.mlst"
		mv "${SAMPDATADIR}/MLST/${isolate_name}.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Pasteur.mlst"
		#Check for "-", unidentified type
		type1=$(tail -n1 ${SAMPDATADIR}/MLST/${isolate_name}_Oxford.mlst | cut -d' ' -f3)
		type2=$(head -n1 ${SAMPDATADIR}/MLST/${isolate_name}_Pasteur.mlst | cut -d' ' -f3)
		if [[ "${type1}" = "-" ]]; then
			singularity -s exec ${src}/singularity_images/srst2.simg getmlst.py --species "Acinetobacter baumannii#1" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
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
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/srst2.simg srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}

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
			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Acinetobacter_baumannii#1-${db_name}.mlst" -t srst2
		fi
		if [[ "${type2}" = "-" ]]; then
			singularity -s exec ${src}/singularity_images/srst2.simg getmlst.py --species "Acinetobacter baumannii#2" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
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
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/srst2.simg srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}

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
			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Acinetobacter_baumannii#2-${db_name}.mlst" -t srst2
		fi
	elif [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
		# Verify that ecoli_2 is default and change accordingly
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/mlst:2.16--0 mlst --scheme "ecoli_2" /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta > "${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst"
		python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst" -t standard
		mv "${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Pasteur.mlst"
		mv "${SAMPDATADIR}/MLST/${isolate_name}.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Achtman.mlst"
		type2=$(tail -n1 ${SAMPDATADIR}/MLST/${isolate_name}_ecoli_2.mlst | cut -d' ' -f3)
		type1=$(head -n1 ${SAMPDATADIR}/MLST/${isolate_name}.mlst | cut -d' ' -f3)
		if [[ "${type1}" = "-" ]]; then
			singularity -s exec ${src}/singularity_images/srst2.simg getmlst.py --species "Escherichia coli#1" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
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
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/srst2.simg srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDir/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}
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

			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Escherichia_coli#1-${db_name}.mlst" -t srst2
		fi
		if [[ "${type2}" = "-" ]]; then
			singularity -s exec ${src}/singularity_images/srst2.simg getmlst.py --species "Escherichia coli#2" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
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
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/srst2.simg srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}
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

			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_Escherichia_coli#2-${db_name}.mlst" -t srst2
		fi
	else
		if [[ "${type}" == "-" ]]; then
			singularity -s exec ${src}/singularity_images/srst2.simg getmlst.py --species "Escherichia coli#1" > "${SAMPDATADIR}/MLST/srst2/getmlst.out"
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
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/srst2.simg srst2 --input_pe /SAMPDIR/srst2/${isolate_name}_S1_L001_R1_001.fastq.gz /SAMPDIR/srst2/${isolate_name}_S1_L001_R2_001.fastq.gz --output /SAMPDIR/srst2/MLST/srst2/${isolate_name} --mlst_db /SAMPDIR/srst2/${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}
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
			python3 "${src}/check_and_fix_MLST.py" -i "${SAMPDATADIR}/MLST/${isolate_name}_srst2_${genus}_${species}-${db_name}.mlst" -t srst2
		fi
		mv "${SAMPDATADIR}/MLST/${isolate_name}.mlst" "${SAMPDATADIR}/MLST/${isolate_name}_Pasteur.mlst"
	fi
	end=$SECONDS
	timeMLST=$((end - start))
	echo "MLST - ${timeMLST} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeMLST))

	exit

	# Try to find any plasmids
	echo "----- Identifying plasmids using plasmidFinder -----"
	start=$SECONDS

	echo "${family}-${genus}"
	# If family is enterobacteriaceae, then run against that DB
	if [[ "${family,}" == "enterobacteriaceae" ]]; then
		echo "Checking against Enterobacteriaceae plasmids"
		#plasmidfinder -i ${SAMPDATADIR}/${inpath} -o ${SAMPDATADIR} -k ${plasmidFinder_identity} -p enterobacteriaceae
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasmidFinder:2.1--0 -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmidFinder -k ${plasmidFinder_identity} -p enterobacteriaceae
		# Rename all files to include ID
		mv ${SAMPDATADIR}/Hit_in_genome_seq.fsa ${SAMPDATADIR}/${isolate_name}_Hit_in_genome_seq_entero.fsa
		mv ${SAMPDATADIR}/Plasmid_seq.fsa ${SAMPDATADIR}/${isolate_name}_Plasmid_seq_enetero.fsa
		mv ${SAMPDATADIR}/results.txt ${SAMPDATADIR}/${isolate_name}_results_entero.txt
		mv ${SAMPDATADIR}/results_tab.txt ${SAMPDATADIR}/${isolate_name}_results_tab_entero.txt
		mv ${SAMPDATADIR}/results_table.txt ${SAMPDATADIR}/${isolate_name}_results_table_summary.txt
	# If family is staph, strp, or enterococcus, then run against the gram positive database
elif [[ "${genus,}" == "staphylococcus" ]] || [[ "${genus,}" == "streptococcus" ]] || [[ "${genus,}" == "enterococcus" ]]; then
		echo "Checking against Staph, Strep, and Enterococcus plasmids"
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasmidFinder:2.1--0 -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmidFinder -k ${plasmidFinder_identity} -p gram_positive
		# Rename all files to include ID
		mv ${SAMPDATADIR}/Hit_in_genome_seq.fsa ${SAMPDATADIR}/${isolate_name}_Hit_in_genome_seq_gramp.fsa
		mv ${SAMPDATADIR}/Plasmid_seq.fsa ${SAMPDATADIR}/${isolate_name}_Plasmid_seq_gramp.fsa
		mv ${SAMPDATADIR}/results.txt ${SAMPDATADIR}/${isolate_name}_results_gramp.txt
		mv ${SAMPDATADIR}/results_tab.txt ${SAMPDATADIR}/${isolate_name}_results_tab_gramp.txt
		mv ${SAMPDATADIR}/results_table.txt ${SAMPDATADIR}/${isolate_name}_results_table_summary.txt
	# Family is not one that has been designated by the creators of plasmidFinder to work well, but still attempting to run against both databases
	else
		echo "Checking against ALL plasmids, but unlikely to find anything"
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasmidFinder:2.1--0 -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmidFinder -k ${plasmidFinder_identity} -p enterobacteriaceae
		# Rename all files to include ID
		mv ${SAMPDATADIR}/Hit_in_genome_seq.fsa ${SAMPDATADIR}/${isolate_name}_Hit_in_genome_seq_entero.fsa
		mv ${SAMPDATADIR}/Plasmid_seq.fsa ${SAMPDATADIR}/${isolate_name}_Plasmid_seq_enetero.fsa
		mv ${SAMPDATADIR}/results.txt ${SAMPDATADIR}/${isolate_name}_results_entero.txt
		mv ${SAMPDATADIR}/results_tab.txt ${SAMPDATADIR}/${isolate_name}_results_tab_entero.txt
		mv ${SAMPDATADIR}/results_table.txt ${SAMPDATADIR}/${isolate_name}_results_table_entero.txt
		#plasmidfinder -i ${SAMPDATADIR}/${inpath} -o ${SAMPDATADIR} -k ${plasmidFinder_identity} -p gram_positive
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasmidFinder:2.1--0 -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmidFinder -k ${plasmidFinder_identity} -p gram_positive
		# Rename all files to include ID
		mv ${SAMPDATADIR}/Hit_in_genome_seq.fsa ${SAMPDATADIR}/${isolate_name}_Hit_in_genome_seq_gramp.fsa
		mv ${SAMPDATADIR}/Plasmid_seq.fsa ${SAMPDATADIR}/${isolate_name}_Plasmid_seq_gramp.fsa
		mv ${SAMPDATADIR}/results.txt ${SAMPDATADIR}/${isolate_name}_results_gramp.txt
		mv ${SAMPDATADIR}/results_tab.txt ${SAMPDATADIR}/${isolate_name}_results_tab_gramp.txt
		mv ${SAMPDATADIR}/results_table.txt ${SAMPDATADIR}/${isolate_name}_results_table_gramp.txt
		cat	${SAMPDATADIR}/${isolate_name}_results_table_gramp.txt ${SAMPDATADIR}/${isolate_name}_results_table_entero.txt > ${SAMPDATADIR}/${isolate_name}_results_table_summary.txt
	fi
	end=$SECONDS
	timeplasfin=$((end - start))
	echo "plasmidFinder - ${timeplasfin} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeplasfin))

	# Run plasFlow if isolate is from the Enterobacteriaceae family  ##### When should we check if this will be expanded?
	if [[ "${family}" == "Enterobacteriaceae" ]]; then
		start=$SECONDS
		${src}/run_plasFlow.sh "${isolate_name}" "${PROJECT}"
		# Create output directory


		if [[ ! -d "${SAMPDATADIR}/plasFlow" ]]; then
			mkdir "${SAMPDATADIR}/plasFlow"
		fi
		if [[ -s "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta" ]]; then
			# Trim contigs a little to 2000 and larger and put through plasflow.
			python3 "${src}/removeShortContigs.py" -i "${SAMPDATADIR}/Assembly/scaffolds.fasta" -t 2000 -s "normal_SPAdes"
			# Renames headers of fasta files accordingly
			python3 "${src}/fasta_headers.py" -i "${SAMPDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta" -o "${SAMPDATADIR}/plasFlow/${isolate_name}_scaffolds_trimmed_2000.fasta"
			# Removes intermeidate fasta file
			rm -r "${SAMPDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta"
			# Run plasflow on newly trimmed assembly file
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasflow:1.1.0--py35_0 --input /SAMPDIR/plasFlow/${isolate_name}_scaffolds_trimmed_2000.fasta --output /SAMPDIR/plasFlow/${isolate_name}_plasFlow_results.tsv --threshold 0.7
			mkdir ${SAMPDATADIR}/plasFlow/bowtie2-index/
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/bowtie2-2.2.9-biocontainers.simg bowtie2-build -f /SAMPDIR/plasFlow/${isolate_name}_plasFlow_results.tsv_chromosomes.fasta /SAMPDIR/plasFlow/bowtie2-index/bowtie2_${isolate_name}_chr
			mkdir ${SAMPDATADIR}/plasFlow/filtered_reads_70/
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/bowtie2-2.2.9-biocontainers.simg bowtie2 -x /SAMPDIR/plasFlow/bowtie2-index/bowtie2_${isolate_name}_chr -1 /SAMPDIR/trimmed/${isolate_name}_R1_001.paired.fq -2 /SAMPDIR/trimmed/${isolate_name}_R2_001.paired.fq -S /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}.sam -p ${procs} --local
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/samtools:1.10--h9402c20_2 samtools view -bS /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}.sam > "${SAMPDATADIR}/plasFlow/filtered_reads_70/${isolate_name}.bam"
  		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/bedtools:2.29.2--hc088bd4_0 bam2fastq --no-aligned -o /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}_R#_bacterial.fastq /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}.bam
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/unicycler:0.4.4--py37h8b12597_2 unicycler -1 /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}_R_1_bacterial.fastq -2 /SAMPDIR/plasFlow/filtered_reads_70/${isolate_name}_R_2_bacterial.fastq -o /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly
			mv "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/assembly.fasta" "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_original.fasta"
			mv "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/assembly.gfa" "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_assembly.gfa"
			python3 ${src}/fasta_headers_plasFlow.py -i "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_original.fasta" -o "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly.fasta"
			python3 "${src}/removeShortContigs.py" -i "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly.fasta" -t 500 -s "plasFlow"
			mv "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly.fasta.TRIMMED.fasta" "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta"
			rm "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly.fasta"
		else
			echo "${SAMPDATADIR}/Assembly/${isolate_name}_scaffolds_trimmed.fasta (Assembly) not found, cant do anything"
		fi


		if [[ -s "${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta" ]]; then
			if [ ! -d "${SAMPDATADIR}/Assembly_Stats_plasFlow" ]; then
				echo "Creating ${SAMPDATADIR}/Assembly_Stats_plasFlow"
		 		mkdir -p "${SAMPDATADIR}/Assembly_Stats_plasFlow"
		 	fi
		 	singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/QUAST5.simg python3 /quast/quast.py -o /SAMPDIR/Assembly_Stats_plasFlow /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta

		 	mv "${SAMPDATADIR}/Assembly_Stats_plasFlow/report.txt" "${SAMPDATADIR}/Assembly_Stats_plasFlow/${isolate_name}_report.txt"
		 	mv "${SAMPDATADIR}/Assembly_Stats_plasFlow/report.tsv" "${SAMPDATADIR}/Assembly_Stats_plasFlow/${isolate_name}_report.tsv"
		 else
			echo "No plasFlow assembly (${SAMPDATADIR}/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta)"
		fi

		# Creates the output c-sstar folder if it does not exist yet
		if [ ! -d "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}" ]; then  #create outdir if absent
			echo "Creating ${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}"
			mkdir -p "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}"
		fi
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/cSSTAR.simg python3 /cSSTAR/c-SSTAR_${csstar_gapping}.py -g /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta -s "${cpsim}" -d "${ResGANNCBI_srst2}" > "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar"

		###################################### FIND WAY TO CATCH FAILURE !!!!!!!!!! ###############################

		# Goes through ResGANNCBI outfile and adds labels as well as resistance conferred to the beginning of the line
		# Takes .sstar file in and outputs as .sstar_grouped
		while IFS= read -r line || [ -n "$line" ]; do

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
			# Catch instances where match length is longer than gene (gaps cause extension)
			#if [[ ${len1} -ge ${len2} ]]; then
			#	plen=100
			# Determine % length match
			#else
			#	plen=$( echo "${len1} ${len2}" | awk '{ printf "%d", ($1*100)/$2 }' )
			#fi
			# Check and display any flags found, otherwise mark it as normal
			if [[ -z "${info1}" ]]; then
				info1="normal"
			else
				info1=${info1::-1}
			fi
			#printf "%-10s %-50s %-15s %-25s %-25s %-40s %-4s %-5d %-5d %-5d\\n" "${source}1" "${resistance}2" "${label1}3" "${info1}4" "${label2}5" "${contig}A" "${percent}B" "${len1}C" "${len2}D" "${plen}E"
			echo "${source}	${resistance}	${label1}	${info1}	${label2}	${contig}	${percent}	${len1}	${len2}	${plen}"
		done < "${SAMPDATADIR}/csstar-plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar" > "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar_grouped"
		# Writes all AR genes to file based on %ID, %length, and finally length of gene
		sort -k7,7nr -k10,10nr -k8,8n "${SAMPDATADIR}/c-sstar_plasFlow/${ResGANNCBI_srst2_filename}_${csstar_gapping}/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}.sstar_grouped" > "${SAMPDATADIR}/c-sstar_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}_sstar_summary.txt"

		# Catches an empty or missing file, adding that no AMR genes were found if no file was created
		if [ ! -s "${SAMPDATADIR}/c-sstar_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}_sstar_summary.txt" ]; then
			echo "No anti-microbial genes were found using c-SSTAR with both resFinder and ARG-ANNOT DBs" > "${SAMPDATADIR}/c-sstar_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.${csstar_gapping}_${cpsim}_sstar_summary.txt"
		fi

		# If family is enterobacteriaceae, then run against that DB
		if [[ "${family,}" == "enterobacteriaceae" ]]; then
			echo "Checking against Enterobacteriaceae plasmids"
			#plasmidfinder -i ${SAMPDATADIR}/${inpath} -o ${SAMPDATADIR} -k ${plasmidFinder_identity} -p enterobacteriaceae
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasmidFinder:2.1--0 -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmid_on_plasFlow -k ${plasmidFinder_identity} -p enterobacteriaceae
			# Rename all files to include ID
			mv ${SAMPDATADIR}/Hit_in_genome_seq.fsa ${SAMPDATADIR}/${isolate_name}_Hit_in_genome_seq_entero.fsa
			mv ${SAMPDATADIR}/Plasmid_seq.fsa ${SAMPDATADIR}/${isolate_name}_Plasmid_seq_enetero.fsa
			mv ${SAMPDATADIR}/results.txt ${SAMPDATADIR}/${isolate_name}_results_entero.txt
			mv ${SAMPDATADIR}/results_tab.txt ${SAMPDATADIR}/${isolate_name}_results_tab_entero.txt
			mv ${SAMPDATADIR}/results_table.txt ${SAMPDATADIR}/${isolate_name}_results_table_summary.txt
		# If family is staph, strp, or enterococcus, then run against the gram positive database
	elif [[ "${genus,}" == "staphylococcus" ]] || [[ "${genus,}" == "streptococcus" ]] || [[ "${genus,}" == "enterococcus" ]]; then
			echo "Checking against Staph, Strep, and Enterococcus plasmids"
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasmidFinder:2.1--0 -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmid_on_plasFlow -k ${plasmidFinder_identity} -p gram_positive
			# Rename all files to include ID
			mv ${SAMPDATADIR}/Hit_in_genome_seq.fsa ${SAMPDATADIR}/${isolate_name}_Hit_in_genome_seq_gramp.fsa
			mv ${SAMPDATADIR}/Plasmid_seq.fsa ${SAMPDATADIR}/${isolate_name}_Plasmid_seq_gramp.fsa
			mv ${SAMPDATADIR}/results.txt ${SAMPDATADIR}/${isolate_name}_results_gramp.txt
			mv ${SAMPDATADIR}/results_tab.txt ${SAMPDATADIR}/${isolate_name}_results_tab_gramp.txt
			mv ${SAMPDATADIR}/results_table.txt ${SAMPDATADIR}/${isolate_name}_results_table_summary.txt
		# Family is not one that has been designated by the creators of plasmidFinder to work well, but still attempting to run against both databases
		else
			echo "Checking against ALL plasmids, but unlikely to find anything"
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasmidFinder:2.1--0 -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmid_on_plasFlow -k ${plasmidFinder_identity} -p enterobacteriaceae
			# Rename all files to include ID
			mv ${SAMPDATADIR}/Hit_in_genome_seq.fsa ${SAMPDATADIR}/${isolate_name}_Hit_in_genome_seq_entero.fsa
			mv ${SAMPDATADIR}/Plasmid_seq.fsa ${SAMPDATADIR}/${isolate_name}_Plasmid_seq_enetero.fsa
			mv ${SAMPDATADIR}/results.txt ${SAMPDATADIR}/${isolate_name}_results_entero.txt
			mv ${SAMPDATADIR}/results_tab.txt ${SAMPDATADIR}/${isolate_name}_results_tab_entero.txt
			mv ${SAMPDATADIR}/results_table.txt ${SAMPDATADIR}/${isolate_name}_results_table_entero.txt
			#plasmidfinder -i ${SAMPDATADIR}/${inpath} -o ${SAMPDATADIR} -k ${plasmidFinder_identity} -p gram_positive
			singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR docker://quay.io/biocontainers/plasmidFinder:2.1--0 -i /SAMPDIR/Assembly/${isolate_name}_scaffolds_trimmed.fasta -o /SAMPDIR/plasmid_on_plasFlow -k ${plasmidFinder_identity} -p gram_positive
			# Rename all files to include ID
			mv ${SAMPDATADIR}/Hit_in_genome_seq.fsa ${SAMPDATADIR}/${isolate_name}_Hit_in_genome_seq_gramp.fsa
			mv ${SAMPDATADIR}/Plasmid_seq.fsa ${SAMPDATADIR}/${isolate_name}_Plasmid_seq_gramp.fsa
			mv ${SAMPDATADIR}/results.txt ${SAMPDATADIR}/${isolate_name}_results_gramp.txt
			mv ${SAMPDATADIR}/results_tab.txt ${SAMPDATADIR}/${isolate_name}_results_tab_gramp.txt
			mv ${SAMPDATADIR}/results_table.txt ${SAMPDATADIR}/${isolate_name}_results_table_gramp.txt
			cat	${SAMPDATADIR}/${isolate_name}_results_table_gramp.txt ${SAMPDATADIR}/${isolate_name}_results_table_entero.txt > ${SAMPDATADIR}/${isolate_name}_results_table_summary.txt
		fi

		if [ ! -d "${SAMPDATADIR}/GAMA_plasFlow" ]; then  #create outdir if absent
			echo "Creating ${SAMPDATADIR}/GAMA_plasFlow"
			mkdir -p "${SAMPDATADIR}/GAMA_plasFlow"
		fi
		singularity -s exec -B ${SAMPDATADIR}:/SAMPDIR ${src}/singularity_images/GAMA_quaisar.simg python3 /GAMA/GAMA_quaisar.py -i /SAMPDIR/plasFlow/Unicycler_assemblies/${isolate_name}_uni_assembly/${isolate_name}_plasmid_assembly_trimmed.fasta -d "${ResGANNCBI_srst2}" -o /SAMPDIR/GAMA_plasFlow/${isolate_name}.${ResGANNCBI_srst2_filename}.GAMA


		end=$SECONDS
		timeplasflow=$((end - start))
		echo "plasmidFlow - ${timeplasflow} seconds" >> "${time_summary_redo}"
		totaltime=$((totaltime + timeplasflow))
	fi

	"${src}/validate_piperun.sh" "${isolate_name}" "${PROJECT}" > "${SAMPDATADIR}/${isolate_name}_pipeline_stats.txt"
	"${src}/sample_cleaner.sh" "${isolate_name}" "${PROJECT}"


	# Extra dump cleanse in case anything else failed
	if [ -n "$(find "${src}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
		echo "Found core dump files at end of processing ${isolate_name} and attempting to delete"
		find "${src}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
	fi

	global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

	# Append total time to bottom of time summary
	echo "Total time: ${totaltime} seconds" >> "${time_summary}"
	echo "Completed at ${global_end_time}"

	# Designate end of this sample #
	echo "

					End of sample ${isolate_name}
					completed at ${global_end_time}

	"

done

# Concatenates lists if this run was an addition to an already processed folder
if [[ -f "${PROJDATADIR}/${PROJECT}_list_original.txt" ]]; then
	cat "${PROJDATADIR}/${PROJECT}_list.txt" >> "${PROJDATADIR}/${PROJECT}_list_original.txt"
	rm "${PROJDATADIR}/${PROJECT}_list.txt"
	mv "${PROJDATADIR}/${PROJECT}_list_original.txt" "${PROJDATADIR}/${PROJECT}_list.txt"
fi

# Run the Seqlog creator on the proper file
declare -A mmb_bugs
while IFS= read -r bug_lines  || [ -n "$bug_lines" ]; do
	bug_genus=$(echo "${bug_lines}" | cut -d'	' -f1)
	bug_species=$(echo "${bug_lines}" | cut -d'	' -f2)
	bug_info=$(echo "${bug_lines}" | cut -d'	' -f3-)
	bug_size=$(echo "${bug_lines}" | cut -d'	' -f6)
	bug_name="${bug_genus:0:1}.${bug_species}"
	#echo "Should be adding ${bug_size} for ${bug_name}"
	mmb_bugs["${bug_name}"]="${bug_size}"
done < ${local_DBs}/MMB_Bugs.txt

# Set output folder as directory of input list
> "${PROJDATADIR}/Seqlog_output.txt"

# Goes through each item on the list and pulls all relevant info
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
		while IFS= read -r line  || [ -n "$line" ]; do
			first=${line::1}
			if [ "${first}" = "S" ]
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
		while IFS= read -r line  || [ -n "$line" ]; do
			first=${line::1}
			if [ "${first}" = "S" ]
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
		line=$(head -n 1 "${SAMPDATADIR}/preQCcounts/${isolate_name}_counts.txt")
		IFS='	' read -r -a qcs <<< "${line}"
		read_qc_info=${qcs[@]:1}
	fi

	source_call=$(head -n1 "${SAMPDATADIR}/${isolate_name}.tax")
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line:0:1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			dec_species=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "G" ]
		then
			dec_genus=$(echo "${line}" | awk -F ' ' '{print $2}')
		fi
	done < "${SAMPDATADIR}/${isolate_name}.tax"

	# Pulls contig info from toms qc analysis file
	contig_info="0(0)\\t0\tNot_in_DB"
		if [[ -s "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv" ]]; then
		counter=0
		while IFS= read -r line  || [ -n "$line" ]; do
			if [ ${counter} -eq 0 ]
			then
				num_contigs=$(sed -n '14p' "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv"| sed -r 's/[\t]+/ /g' | cut -d' ' -f3 )
			elif [ ${counter} -eq 1 ]
			then
				ass_length=$(sed -n '16p' "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
				#Check Assembly ratio against expected size to see if it is missing a large portion or if there is contamination/double genome
				dec_genus_initial="${dec_genus:0:1}"
				if [[ "${dec_genus_initial}" = "[" ]] || [[ "${dec_genus_initial}" = "(" ]]; then
					dec_genus_initial="${dec_genus:1:1}"
				fi
				ass_ID="${dec_genus_initial}.${dec_species}"
				echo "About to check mmb_bugs[${ass_ID}]"
				if [[ ! -z "${mmb_bugs[${ass_ID}]}" ]]; then
					#echo "Found Bug in DB: ${ass_ID}-${mmb_bugs[${ass_ID}]}"
					ass_ratio=$(awk -v p="${ass_length}" -v q="${mmb_bugs[${ass_ID}]}" 'BEGIN{printf("%.2f",p/q)}')
				else
					ass_ratio="Not_in_DB"
				fi
			elif [ ${counter} -eq 3 ]
			then
				N50=$(sed -n '18p' "${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv"  | sed -r 's/[\t]+/ /g'| cut -d' ' -f2)
			fi
			counter=$((counter+1))
		done < ${SAMPDATADIR}/Assembly_Stats/${isolate_name}_report.tsv
		contig_info=$(echo -e "${num_contigs}\\t${ass_length}\\t${ass_ratio}")
	fi

	# Pulls busco info from summary file
	busco_info="No BUSCO performed"
	if [[ -s "${SAMPDATADIR}/BUSCO/short_summary_${isolate_name}.txt" ]]; then
		while IFS= read -r line  || [ -n "$line" ]; do
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
	# Rename files in old formating convention
	if [[ -s "${SAMPDATADIR}/ANI/best_hits_ordered.txt" ]]; then
		mv "${SAMPDATADIR}/ANI/best_hits_ordered.txt" "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_${genus}).txt"
	fi
	# If 1 and only 1 file exists pull the first line as the best hit information
	if [[ -s "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_All.txt" ]]; then
		ani_info=$(head -n 1 "${SAMPDATADIR}/ANI/best_ANI_hits_ordered(${isolate_name}_vs_All).txt")
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
	if [[ ${ass_length} -gt 0 ]]; then
		avg_coverage=$(bc <<<"scale=2 ; ${q30_reads} / ${ass_length}")
	else
		avg_coverage="N/A"
	fi
	# Replace all spaces in qc_info as tabs
	read_qc_info=$(echo "${read_qc_info}" | tr '[:blank:]' \\t)

	# Get the date to show when the log was made
	NOW=$(date +"%m/%d/%Y")

	# Add all pertinent info to the output file in the correct formatting to add to MMB_Seq log
	echo -e "${isolate_name}\\t${NOW}\\t${g_s_reads}\\t${g_s_assembled}\\t${g_s_16s}\\t${read_qc_info}\\t${avg_coverage}\\t${contig_info}\\t${busco_info}\\t${ani_info}\\r" >> "${output_folder}/Seqlog_output.txt"
done < ${list_path}

# Get run summary info to send in an email
runsumdate=$(date "+%m_%d_%Y_at_%Hh_%Mm")
${src}/run_sum.sh ${PROJECT}

# Copy the config file to the log directory so as not to hold up any future quaisar runs that count the number of config files present, but for some reason does not remove from script folder
echo "Moving config file(${config_file}) to log directory ($log_dir)"
	mv "${config_file}" "${log_dir}/config_${PROJECT}.sh"

end_date=$(date "+%m_%d_%Y_at_%Hh_%Mm")
echo "Run ended at ${end_date}" >> "${log_dir}/${PROJECT}_on_${run_start_time}.log"
