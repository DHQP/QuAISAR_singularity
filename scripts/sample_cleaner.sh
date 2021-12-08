#!/bin/bash -l

#$ -o sample_cleaner.out
#$ -e sample_cleaner.err
#$ -N sample_cleaner
#$ -cwd
#$ -q short.q

#
# Description: Script uses Gulviks SPAdes cleaner along with general folder cleanup to decrease footprint of samples after processing
#
# Usage ./sample_cleaner.sh   path_to_sample_folder
#
# Output location: No output created
#
# Modules required: None
#
# v1.0.1 (05/14/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./sample_cleaner.sh  sample_name MiSeq_Run_ID"
	echo "Will clean ${1} folder"
	exit 0
elif [[ ! -d "${1}" ]]; then
	echo "Sample folder (${1}) does not exist"
	exit 2
else
	echo "Cleaning ${1}"
fi

# Get path of this script so any other scripts calls can be referential
self_path=$( realpath "$0"  ) && dirname "$self_path"

# Set main sample folder to clean
sample_folder="${1}"
echo "Source - ${sample_folder}"
sample_name=$(echo "${sample_folder}" | rev | cut -d'/' -f1 | rev)
echo "Sample_ID - ${sample_name}"
# Remove the localANIDB from ANI output folder, if found
echo "Cleaning ANI"
if [ -d "${sample_folder}/ANI/localANIDB" ]; then
	echo "removing localANIDB"
	rm -r "${sample_folder}/ANI/localANIDB"
fi
if [ -d "${sample_folder}/ANI/localANIDB_REFSEQ" ]; then
	echo "removing localANIDB_REFSEQ"
	rm -r "${sample_folder}/ANI/localANIDB_REFSEQ"
fi
if [ -d "${sample_folder}/ANI/localANIDB_full" ]; then
	echo "removing localANIDB_full"
	rm -r "${sample_folder}/ANI/localANIDB_full"
fi
if [ -d "${sample_folder}/ANI/temp" ]; then
	echo "removing temp"
	rm -r "${sample_folder}/ANI/temp"
fi
#if [ -d "${sample_folder}/ANI/aniM_REFSEQ" ]; then
#	echo "removing nucmers"
#	rm -r "${sample_folder}/ANI/aniM_REFSEQ/nucmer_output.tar.gz"
#fi
# Remove the hmmer output from the BUSCO folder
echo "Cleaning BUSCO"
if [ -d "${sample_folder}/BUSCO/hmmer_output" ]; then
	echo "removing hmmer output"
	rm -r "${sample_folder}/BUSCO/hmmer_output"
fi
# Use Gulviks cleaner script on regular SPAdes output
echo "Cleaning Assembly Folder"
if [ -d "${sample_folder}/Assembly" ]; then
	echo "Using Gulviks SPAdes cleaner on Assembly"
	${self_path}/gulvic_SPAdes_cleaner.sh "${sample_folder}/Assembly"
	rm "${sample_folder}/Assembly/${sample_name}_scaffolds_trimmed_original.fasta"
fi

# Cleaning Assembly Stats folder of extra unused files made by QUAST
echo "Cleaning Assembly Stats Folder"
if [ -d "${sample_folder}/Assembly_Stats" ]; then
	rm "${sample_folder}/Assembly_Stats/${sample_name}_report.txt"
	rm "${sample_folder}/Assembly_Stats/report.tex"
	rm "${sample_folder}/Assembly_Stats/transposed_report.t"*
	rm "${sample_folder}/Assembly_Stats/icarus.html"
	rm -r "${sample_folder}/Assembly_Stats/icarus_viewers"
	rm -r "${sample_folder}/Assembly_Stats/basic_stats"
fi

# # Cleaning Assembly Stats plasFlow folder of extra unused files made by QUAST
# echo "Cleaning Assembly Stats plasFlow Folder"
# if [ -d "${sample_folder}/Assembly_Stats_plasFlow" ]; then
# 	rm "${sample_folder}/Assembly_Stats_plasFlow/${sample_name}_report.txt"
# 	rm "${sample_folder}/Assembly_Stats_plasFlow/report.tex"
# 	rm "${sample_folder}/Assembly_Stats_plasFlow/transposed_report.t"*
# 	rm "${sample_folder}/Assembly_Stats_plasFlow/icarus.html"
# 	rm -r "${sample_folder}/Assembly_Stats_plasFlow/icarus_viewers"
# 	rm -r "${sample_folder}/Assembly_Stats_plasFlow/basic_stats"
# fi

# Clean kraken folder (only for reads as it is MUCH larger)
if [ -d "${sample_folder}/kraken/preAssembly" ]; then
	gzip ${sample_folder}/kraken/preAssembly/${sample_name}_paired.classified
	#rm ${sample_folder}/kraken/preAssembly/${sample_name}_paired.classified
	gzip ${sample_folder}/kraken/preAssembly/${sample_name}_paired.kraken
fi

# Clean kraken folder (only for reads as it is MUCH larger)
if [ -d "${sample_folder}/kraken/postAssembly" ]; then
	gzip ${sample_folder}/kraken/postAssembly/${sample_name}_assembled.classified
	#rm ${sample_folder}/kraken/postAssembly/${sample_name}_assembled.classified
	gzip ${sample_folder}/kraken/postAssembly/${sample_name}_assembled.kraken
fi

# # Clean plasFlow folder of filtered reads
# if [[ -d "${sample_folder}/plasFlow/filtered_reads_70" ]]; then
# 	rm -r "${sample_folder}/plasFlow/filtered_reads_70"
# fi

# Clean plasmidFinder folder
if [[ -d "${sample_folder}/plasmidFinder/tmp" ]]; then
	rm -r "${sample_folder}/plasmidFinder/tmp"
fi

# # Clean plasmidFinder folder
# if [[ -d "${sample_folder}/plasmidFinder_on_plasFlow/tmp" ]]; then
# 	rm -r "${sample_folder}/plasmidFinder_on_plasFlow/tmp"
# fi

# Remove hmm_output folder from BUSCO analysis if found
echo "Cleaning BUSCO Folder"
if [ -d "${sample_folder}/BUSCO/hmm_output" ]; then
	echo "Deleting hmm_output"
	rm -r "${sample_folder}/BUSCO/hmm_output"
fi
# echo "Cleaning GOTTCHA Folder"
# # Remove splitrim fodler from gottcha output, if found
# if [ -d "${sample_folder}/gottcha/gottcha_S/${sample_name}_temp/splitrim" ]; then
# 	echo "Deleting splitrim folder"
# 	rm -r "${sample_folder}/gottcha/gottcha_S/${sample_name}_temp/splitrim"
# fi
# Removed intermediate folder that has reads with no adapters, but have not been trimmed yet
echo "Cleaning Adapter Folder"
if [ -d "${sample_folder}/removedAdapters" ]; then
	echo "Deleting adapterless reads"
	rm -r "${sample_folder}/removedAdapters/"*".fsq"
fi
# Clean trimmed folder of catted and unpaired reads, leaving only R1 and R2
echo "Cleaning Trimmed Folder"
if [ -d "${sample_folder}/trimmed" ]; then
	echo "Deleting extraneous reads"
	if [ -f "${sample_folder}/trimmed/${sample_name}.paired.fq" ]; then
		echo "Deleting catted paired reads"
		rm "${sample_folder}/trimmed/${sample_name}.paired.fq"
	fi
	#if [ -f "${sample_folder}/trimmed/${sample_name}.single.fq" ]; then
	#	echo "Deleting catted single reads"
	#	rm "${sample_folder}/trimmed/${sample_name}.single.fq"
	#fi
	if [ -f "${sample_folder}/trimmed/${sample_name}_R1_001.unpaired.fq" ]; then
		echo "Deleting unpaired R1 reads"
		rm "${sample_folder}/trimmed/${sample_name}_R1_001.unpaired.fq"
	fi
	if [ -f "${sample_folder}/trimmed/${sample_name}_R2_001.unpaired.fq" ]; then
		echo "Deleting unpaired R2 reads"
		rm "${sample_folder}/trimmed/${sample_name}_R2_001.unpaired.fq"
	fi
fi
# Clean FASTQ folder by zipping any unzipped reads
# if [ -s "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ] && [ -s "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz"]; then
# 	echo "Clumping R1 and R2"
#
# 	if [ -f "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq" ]; then
# 		rm "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
# 	fi
# 	if [ -f "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq" ]; then
# 		rm "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
# 	fi
if [ -f "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq" ]; then
	#echo "Found unzipped FASTQ"
	if [ ! -f "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ]; then
		echo "Zipping R1"
		gzip -c "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
		if [ -s "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ]; then
			echo "Zipping seems to have been successful, deleting ${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
			rm -r "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
		fi
	else
		echo "Zipped file found, checking for substance"
		if [ -s "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ]; then
			echo "Zipped file not zero, deleting ${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
			rm -r "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
		else
			echo "Zip file was empty, trying to rezip"
			gzip -c "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
			if [ -s "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ]; then
				echo "Rezip worked,Deleting ${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
				rm -r "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
			else
				echo "Rezip did not work, what to do now???"
			fi
		fi
	fi
fi
if [ -f "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq" ]; then
	#echo "Found unzipped FASTQ"
	if [ ! -f "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz" ]; then
		echo "Zipping R2"
		gzip -c "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
		if [ -s "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz" ]; then
			echo "Zipping seems to have been successful, deleting ${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
			rm -r "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
		fi
	else
		echo "Zipped file found, checking for substance"
		if [ -s "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz" ]; then
			echo "Zipped file not zero, deleting ${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
			rm -r "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
		else
			echo "Zip file was empty, trying to rezip"
			gzip -c "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
			if [ -s "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz" ]; then
				echo "Rezip worked,Deleting ${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
				rm -r "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
			else
				echo "Rezip did not work, what to do now???"
			fi
		fi
	fi
fi
# Zip trimmed R1 read, if not already done
if [ -f "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq" ]; then
	#echo "Found unzipped FASTQ"
	if [ ! -f "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
		echo "Zipping paired R1"
		gzip "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
		if [ -s "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
			echo "Zipping seems to have been successful, deleting ${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
			rm -r "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
		fi
	else
		echo "Zipped file found, checking for substance"
		if [ -s "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
			echo "Zipped file not zero, deleting ${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
			rm -r "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
		else
			echo "Zip file was empty, trying to rezip"
			gzip "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
			if [ -s "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
				echo "Zipping seems to have been successful, deleting ${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
				rm -r "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
			else
				echo "Rezip did not work, what to do now???"
			fi
		fi
	fi
fi
# Zip trimmed R2 read, if not already done
if [ -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq" ]; then
	#echo "Found unzipped FASTQ"
	if [ ! -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
		echo "Zipping paired R2"
		gzip "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
		if [ -s "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
			echo "Zipping seems to have been successful, deleting ${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
			rm -r "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
		fi
	else
		echo "Zipped file found, checking for substance"
		if [ -s "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
			echo "Zipped file not zero, deleting ${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
			rm -r "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
		else
			echo "Zip file was empty, trying to rezip"
			gzip "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
			if [ -s "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
				echo "Zipping seems to have been successful, deleting ${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
				rm -r "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
			else
				echo "Rezip did not work, what to do now???"
			fi
		fi
	fi
fi
# Zip single fq files if present
if [ -f "${sample_folder}/trimmed/${sample_name}.single.fq" ]; then
	if [ ! -f "${sample_folder}/trimmed/${sample_name}.single.fq.gz" ]; then
		gzip "${sample_folder}/trimmed/${sample_name}.single.fq"
	fi
fi

# Try to reduce size of zipped files further by using clumpify
#if [[ -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]] && [[ -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]]; then
#	clumpify in1="${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" in2="${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" out1="${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.clumped.gz" out2="${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.clumped.gz" reorder
#fi

# echo "Sample ${1} should now be clean" >> "${output_dir}/cleaned_sample_list.txt"
