#!/bin/sh -l

#$ -o validate_piperun.out
#$ -e validate_piperun.err
#$ -N validate_piperun
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Checking to see if all standard reported sections of a sample have completed successfully
#
# Usage: ./validate_piprun.sh   sample_name   miseq run id [gapping (gapped|ungapped)] [sim (40|80|95|998|99|100)]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/
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
	echo "Empty sample name supplied to validate_piperun.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./validate_piperun.sh   sample_name   miseq_run_ID [gapping (gapped|ungapped)] [sim (40|80|95|998|99|100)]"
	echo "Output is only printed to screen, Pipe to file if desired"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id supplied to validate_piperun.sh, exiting"
	exit 1
fi

OUTDATADIR=${processed}/${2}/${1}

# Creates and prints header info for the sample being processed
today=$(date)
echo "----------Checking ${2}/${1} for successful completion on ----------"
echo "Sample output folder starts at: " "${processed}/${2}/${1}"
status="SUCCESS"
# Checks to see if the sample has a time summary file associated with it
if [[ -s "${OUTDATADIR}/time_summary.txt" ]]; then
	mv "${OUTDATADIR}/time_summary.txt" "${OUTDATADIR}/${1}_time_summary.txt"
fi
printf "%-20s: %-8s : %s\\n" "Summarized" "SUCCESS" "${today}"
if [[ -s "${OUTDATADIR}/${1}_time_summary.txt" ]]; then
	time=$(tail -1 "${OUTDATADIR}/${1}_time_summary.txt" | cut -d' ' -f3)
	printf "%-20s: %-8s : %s\\n" "Time" "SUCCESS" "${time} seconds"
else
	printf "%-20s: %-8s : %s\\n" "Time" "ALERT" "No time summary file found"
	status="ALERT"
fi
#Checking existence of FASTQ files
if [[ -s "${OUTDATADIR}/FASTQs/${1}_R1_001.fastq" ]] && [[ -s "${OUTDATADIR}/FASTQs/${1}_R2_001.fastq" ]]; then
	printf "%-20s: %-8s : %s\\n" "FASTQs" "SUCCESS" "Unzipped"
	:
elif [[ -s "${OUTDATADIR}/FASTQs/${1}_R1_001.fastq" ]]; then
	printf "%-20s: %-8s : %s\\n" "FASTQs" "WARNING" "Only R1 found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
elif [[ -s "${OUTDATADIR}/FASTQs/${1}_R2_001.fastq" ]]; then
	printf "%-20s: %-8s : %s\\n" "FASTQs" "WARNING" "Only R2 found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
elif [[ -s "${OUTDATADIR}/FASTQs/${1}_R1_001.fastq.gz" ]] && [[ -s "${OUTDATADIR}/FASTQs/${1}_R2_001.fastq.gz" ]]; then
	printf "%-20s: %-8s : %s\\n" "FASTQs" "SUCCESS" "Zipped"
	:
elif [[ -s "${OUTDATADIR}/FASTQs/${1}_R1_001.fastq" ]]; then
	printf "%-20s: %-8s : %s\\n" "FASTQs" "WARNING" "Zipped, but only R1 found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
elif [[ -s "${OUTDATADIR}/FASTQs/${1}_R2_001.fastq" ]]; then
	printf "%-20s: %-8s : %s\\n" "FASTQs" "WARNING" "Zipped, but only R2 found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
else
	printf "%-20s: %-8s : %s\\n" "FASTQs" "FAILED" "No reads found"
	status="FAILED"
fi
#Checking QC counts
if [[ -s "${OUTDATADIR}/preQCcounts/${1}_counts.txt" ]]; then
	reads_pre=$(head -n 1 "${OUTDATADIR}/preQCcounts/${1}_counts.txt" | cut -d'	' -f13)
	pairs_pre=$((reads_pre/2))
	Q30_R1=$(head -n 1 "${OUTDATADIR}/preQCcounts/${1}_counts.txt" | cut -d'	' -f10)
	Q30_R1_rounded=$(echo "${Q30_R1}"  | cut -d'.' -f2)
	Q30_R1_rounded=$(echo "${Q30_R1_rounded::2}")
	Q30_R2=$(head -n 1 "${OUTDATADIR}/preQCcounts/${1}_counts.txt" | cut -d'	' -f11)
	Q30_R2_rounded=$(echo "${Q30_R2}"  | cut -d'.' -f2)
	Q30_R2_rounded=$(echo "${Q30_R2_rounded::2}")
	if [[ "${reads_pre}" -le 1000000 ]]; then
		printf "%-20s: %-8s : %s\\n" "QC counts" "WARNING" "Low individual read count before trimming: ${reads_pre} (${pairs_pre} paired reads)"
		status="WARNING"
	else
		printf "%-20s: %-8s : %s\\n" "QC counts" "SUCCESS" "${reads_pre} individual reads found in sample (${pairs_pre} paired reads)"
		if [[ "${Q30_R1_rounded}" -lt 90 ]]; then
			printf "%-20s: %-8s : %s\\n" "Q30_R1%" "WARNING" "Q30_R1% below 90(${Q30_R1_rounded}%)"
			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
				status="WARNING"
			fi
		else
			printf "%-20s: %-8s : %s\\n" "Q30_R1%" "SUCCESS" "Q30_R1% at ${Q30_R1_rounded}% (Threshold is 90)"
		fi
		if [[ "${Q30_R2_rounded}" -lt 70 ]]; then
			printf "%-20s: %-8s : %s\\n" "Q30_R2%" "WARNING" "Q30_R2% below 70(${Q30_R2_rounded}%)"
			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
				status="WARNING"
			fi
		else
			printf "%-20s: %-8s : %s\\n" "Q30_R2%" "SUCCESS" "Q30_R2% at ${Q30_R2_rounded}% (Threshold is 70)"
		fi
	fi
else
	printf "%-20s: %-8s : %s\\n" "QC counts" "FAILED" "/preQCcounts/${1}_counts.txt not found"
	printf "%-20s: %-8s : %s\\n" "Q30_R1%" "FAILED" "/preQCcounts/${1}_counts.txt not found"
	printf "%-20s: %-8s : %s\\n" "Q30_R2%" "FAILED" "/preQCcounts/${1}_counts.txt not found"
	status="FAILED"
fi

### This folder is now deleted afterwards and therefore is no longer checked
#Checking BBDUK output folder
#if [[ -s "${OUTDATADIR}/removedAdapters/${1}-noPhiX-R1.fsq" ]] && [[ -s "${OUTDATADIR}/removedAdapters/${1}-noPhiX-R2.fsq" ]]; then
#	#printf "%-20s: %-8s : %s\\n" "Adapter Removal" "SUCCESS" "Found"
#	:
#elif [[ -s "${OUTDATADIR}/removedAdapters/${1}-noPhiX-R1.fsq" ]]; then
#	printf "%-20s: %-8s : %s\\n" "Adapter Removal" "WARNING" "Only R1 found"
#	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
#		status="WARNING"
#	fi
#elif [[ -s "${OUTDATADIR}/removedAdapters/${1}-noPhiX-R2.fsq" ]]; then
#	printf "%-20s: %-8s : %s\\n" "Adapter Removal" "WARNING" "Only R2 found"
#	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
#		status="WARNING"
#	fi
#else
#	printf "%-20s: %-8s : %s\\n" "Adapter Removal" "FAILED" "/removedAdapters/${1}-noPhiX-R1.fsq & /removedAdapters/${1}-noPhiX-R2.fsq not found"
#	status="FAILED"
#fi

#Checking Trimmomatic output folder
if [[ -s "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" ]] && [[ -s "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" ]]; then
	printf "%-20s: %-8s : %s\\n" "Trimming" "SUCCESS" "Unzipped"
	:
elif [[ -s "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" ]]; then
	printf "%-20s: %-8s : %s\\n" "Trimming" "WARNING" "Only R1 found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
elif [[ -s "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" ]]; then
	printf "%-20s: %-8s : %s\\n" "Trimming" "WARNING" "Only R2 found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
elif [[ -s "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq.gz" ]] && [[ -s "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq.gz" ]]; then
	printf "%-20s: %-8s : %s\\n" "Trimming" "SUCCESS" "Zipped"
	:
elif [[ -s "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq.gz" ]]; then
	printf "%-20s: %-8s : %s\\n" "Trimming" "WARNING" "Zipped, but only R1 found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
elif [[ -s "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq.gz" ]]; then
	printf "%-20s: %-8s : %s\\n" "Trimming" "WARNING" "Zipped, but only R2 found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
else
	printf "%-20s: %-8s : %s\\n" "Trimming" "FAILED" "/trimmed/${1}_R1_001.paired.fq & /trimmed/${1}_R2_001.paired.fq not found"
	status="FAILED"
fi
#Checking QC counts after trimming
if [[ -s "${OUTDATADIR}/preQCcounts/${1}_trimmed_counts.txt" ]]; then
	reads_post=$(head -n 1 "${OUTDATADIR}/preQCcounts/${1}_trimmed_counts.txt" | cut -d'	' -f13)
	pairs_post=$((reads_post/2))
	loss=$(echo "scale=2; 100*(${reads_pre} - ${reads_post}) / ${reads_pre}" | bc )
	if [[ "${reads_post}" -le 500000 ]]; then
		printf "%-20s: %-8s : %s\\n" "QC count after trim" "WARNING" "Low individual read count after trimming: ${reads_post} (${pairs_post} paired reads)"
		status="WARNING"
	else
		printf "%-20s: %-8s : %s\\n" "QC count after trim" "SUCCESS" "${reads_post} individual reads (${pairs_post} paired reads) after trim. ${loss}% loss"
	fi
else
	printf "%-20s: %-8s : %s\\n" "QC count after trim" "FAILED" "/preQCcounts/${1}_trimmed_counts.txt not found"
	status="FAILED"
fi


#Check kraken on preAssembly
kraken_pre_success=false
if [[ -s "${OUTDATADIR}/kraken/preAssembly/${1}_paired.kraken" ]]; then
	#printf "%-20s: %-8s : %s\\n" "kraken preassembly" "SUCCESS" "Found"
	kraken_pre_success=true
else
	printf "%-20s: %-8s : %s\\n" "kraken preassembly" "FAILED" "/kraken/preAssembly/${1}_paired.kraken not found"
	status="FAILED"
fi

#Check Krona output
if [[ "${kraken_pre_success}" = true ]]; then
	if [[ -s "${OUTDATADIR}/kraken/preAssembly/${1}_paired.krona" ]] && [[ -s "${OUTDATADIR}/kraken/preAssembly/${1}_paired.html" ]]; then
		#printf "%-20s: %-8s : %s\\n" "krona-kraken-preasmb" "SUCCESS" "Found"
		:
	else
		printf "%-20s: %-8s : %s\\n" "krona-kraken-preasmb" "FAILED" "/kraken/preAssembly/${1}_paired.krona & /kraken/preAssembly/${1}_paired.html not found"
		status="FAILED"
	fi
else
	printf "%-20s: %-8s : %s\\n" "krona-kraken-preasmb" "FAILED" "preassembly kraken did not complete successfully"
	status="FAILED"
fi

#Check extraction and unclassified value
if [[ -s "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt" ]]; then
	# Extracts many elements of the summary file to report unclassified and species classified reads and percentages
	unclass=$(head -n 1 "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt" | cut -d' ' -f2)
	#true_unclass=$(head -n 1 "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
	domain=$(sed -n '2p' "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt" | cut -d' ' -f2)
	genuspre=$(sed -n '7p' "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt" | cut -d' ' -f4)
	speciespre=$(sed -n '8p' "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt" | cut -d' ' -f4)
	speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt" | cut -d' ' -f2)
	#true_speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
	# If there are no reads at the domain level, then report no classified reads
	if (( $(echo "${domain} <= 0" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Pre Classify" "FAILED" "There are no classified reads (Did pre assembly kraken fail too?)"
		status="FAILED"
	# If there are classified reads then check to see if percent unclassifed falls above the threshold limit. Report warning if too high or success and stats if below
	else
		if (( $(echo "${unclass} > ${unclass_flag}" | bc -l) )); then
			printf "%-20s: %-8s : %s\\n" "Pre Classify" "WARNING" "unclassified reads comprise ${unclass}% of total"
			if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
				status="WARNING"
			fi
		else
			#printf "%-20s: %-8s : %s\\n" "Pre Classify" "SUCCESS" "${speciespercent}%${true_speciespercent%} ${genuspre} ${speciespre} with ${unclass}%${true_unclass%} unclassified reads"
			printf "%-20s: %-8s : %s\\n" "Pre Classify" "SUCCESS" "${speciespercent}% ${genuspre} ${speciespre} with ${unclass}% unclassified reads"
		fi
	fi
# If no summary file was found
else
	printf "%-20s: %-8s : %s\\n" "Pre Classify" "FAILED" "${OUTDATADIR}/kraken/preAssembly/${1}_kraken_summary_paired.txt not found"
	status="FAILED"
fi

# Quick separate check for contamination by finding # of species above ${contamination_threshold} in list file from kraken
if [[ -s "${OUTDATADIR}/kraken/preAssembly/${1}_paired.list" ]]; then
	number_of_species=0
	while IFS= read -r line; do
		arrLine=(${line})
		# First element in array is the percent of reads identified as the current taxa
		percent=${arrLine[0]}
		percent_integer=$(echo "${percent}" | cut -d'.' -f1)
		# 3rd element is the taxon level classification
		# echo "${percent_integer} vs ${contamination_threshold}"
		classification=${arrLine[3]}
		if [[ "${classification}" == "S" ]] && (( percent_integer > contamination_threshold )); then
			#echo "Adding ${arrLine[5]}-${percent_integer}-${contamination_threshold} to list"
			number_of_species=$(( number_of_species + 1 ))
		fi
	done < ${OUTDATADIR}/kraken/preAssembly/${1}_paired.list
	if [[ "${number_of_species}" -gt 1 ]]; then
		printf "%-20s: %-8s : %s\\n" "pre Class Contam." "WARNING" "${number_of_species} species have been found above the ${contamination_threshold}% threshold"
		if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
				status="WARNING"
		fi
	elif [[ "${number_of_species}" -eq 1 ]]; then
		:
	else
		printf "%-20s: %-8s : %s\\n" "pre Class Contam." "FAILED" "No species have been found above the ${contamination_threshold}% threshold"
	fi
	#echo "Number of species: ${number_of_species}"
fi

#Check gottcha_S output for TSV ouput and the krona file
if [[ -s "${OUTDATADIR}/gottcha/gottcha_S/${1}.gottcha_full.tsv" ]] && [[ -s "${OUTDATADIR}/gottcha/${1}_species.krona.html" ]]; then
	#printf "%-20s: %-8s : %s\\n" "GOTTCHA_S" "SUCCESS" "Found"
	:
elif [[ -s "${OUTDATADIR}/gottcha/gottcha_S/${1}.gottcha_full.tsv" ]]; then
	printf "%-20s: %-8s : %s\\n" "GOTTCHA_S" "WARNING" "No Krona output found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
elif [[ -s "${OUTDATADIR}/gottcha/${1}_species.krona.html" ]]; then
	printf "%-20s: %-8s : %s\\n" "GOTTCHA_S" "WARNING" "No TSV file found"
	if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		status="WARNING"
	fi
else
	printf "%-20s: %-8s : %s\\n" "GOTTCHA_S" "FAILED" "/gottcha/gottcha_S/${1}.gottcha_full.tsv & /gottcha/${1}_species.krona.html not found"
	status="FAILED"
fi

#Check extraction of gottcha id
if [[ -s "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt" ]]; then
	# Extracts many elements of the summary file to report unclassified and species classified reads and percentages
	unclass=$(head -n 1 "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt" | cut -d' ' -f2)
	#true_unclass=$(head -n 1 "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt" | cut -d' ' -f3) # | sed -r 's/[)]+/%)/g')
	phylumpercent=$(sed -n '3p' "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt" | cut -d' ' -f2)
	genuspre=$(sed -n '7p' "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt"| cut -d' ' -f4)
	speciespre=$(sed -n '8p' "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt" | cut -d' ' -f5)
	speciespercent=$(sed -n '8p' "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt" | cut -d' ' -f2)
	true_speciespercent=$(sed -n '8p' "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
	# Gottcha only classifies up to phylum and therefore if no phylum reads, there are no reads
	if (( $(echo "${phylumpercent} <= 0" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "GottchaV1 Classifier" "FAILED" "There are no classified reads"
		status="FAILED"
	# If there are phylum level reads then check to see the percentage. If it falls below the threshold (set in config.sh) report it as a warning, otherwise report all necessary stats
	else
		if (( $(echo "${unclass} > ${unclass_flag}" | bc -l) )); then
			printf "%-20s: %-8s : %s\\n" "GottchaV1 Classifier" "WARNING" "unclassified reads comprise ${unclass}% of total"
			if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
				status="WARNING"
			fi
		else
			printf "%-20s: %-8s : %s\\n" "GottchaV1 Classifier" "SUCCESS" "${speciespercent}%${true_speciespercent} ${genuspre} ${speciespre} with ${unclass}% unclassified reads"
		fi
	fi
# If the summary file does not exist, report as such
else
	printf "%-20s: %-8s : %s\\n" "GottchaV1 Classifier" "FAILED" "${OUTDATADIR}/gottcha/${1}_gottcha_species_summary.txt not found"
	status="FAILED"
fi

# Quick separate check for contamination by finding # of species above ${contamination_threshold} in list file from kraken
if [[ -s "${OUTDATADIR}/gottcha/gottcha_S/${1}.gottcha.tsv" ]]; then
	number_of_species=0
	while IFS= read -r line; do
		# Convert the perfect match to proper format from 1.00 to 100
		if [[ "${line[2]}" = "1.0000" ]] || [[ "${line[2]}" -eq 1 ]]; then
			percent_integer=100
		# Convert all non-perfect matches to the correct matching percent values
		else
			percent="${line[2]:2:2}.${line[2]:4:2}"
			percent_integer=$(echo "${percent}" | cut -d'.' -f1)
		fi
		# Convert a no-match to the correct percent value
		if [[ "${percent}" = "00.00" ]]; then
			percent_integer=0
		fi
		# Takes the first letter of the first column as shorthand for identifying the taxonomic level
		classification="${line[0]::1}"
		if [[ "${classification}" == "s" ]] && (( percent_integer > contamination_threshold )); then
			number_of_species=$(( number_of_species + 1 ))
		fi
	done < ${OUTDATADIR}/gottcha/gottcha_S/${1}.gottcha.tsv
	if [[ $number_of_species -gt 1 ]]; then
		# Holding off on putting a cutoff here, as we cant tell what is an acceptable value to use
		#printf "%-20s: %-8s : %s\\n" "gottcha Contam." "WARNING" "${number_of_species} species have been found above the ${contamination_threshold}% threshold"
		#if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
		#	status="WARNING"
		#fi
		:
	elif [[ "${number_of_species}" -eq 1 ]]; then
		:
	else
		# Holding off on putting a cutoff here, as we cant tell what is an acceptable value to use
		#printf "%-20s: %-8s : %s\\n" "gottcha Contam." "FAILED" "No species have been found above the ${contamination_threshold}% threshold"
		:
	fi
	#echo "Number of species: ${number_of_species}"
fi

#Check spades assembly
if [[ -s "${OUTDATADIR}/Assembly/scaffolds.fasta" ]]; then
	# Count the number of '>' in the assembly file before trimming
	full_scaffolds=">"
	full_scaffolds=$(grep -c ${full_scaffolds} "${OUTDATADIR}/Assembly/scaffolds.fasta")
	printf "%-20s: %-8s : %s\\n" "Assembly" "SUCCESS" "${full_scaffolds} scaffolds found"
else
	printf "%-20s: %-8s : %s\\n" "Assembly" "FAILED" "/Assembly/scaffolds.fasta not found"
	status="FAILED"
fi

# #Check plasFlow plasmid assembly
plasmidsFoundviaplasFlow=0
if [[ -d "${OUTDATADIR}/plasFlow" ]]; then
	if [[ -s "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_original.fasta" ]]; then
		# Count the number of '>' in the assembly file before trimming
		plas_scaffolds=">"
		plas_scaffolds=$(grep -c ${plas_scaffolds} "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_original.fasta")
		if [ -z ${plas_scaffolds} ]; then
			plas_scaffolds=0
		fi
		if [[ "${plas_scaffolds}" -gt 0 ]]; then
			printf "%-20s: %-8s : %s\\n" "plasmid Assembly" "SUCCESS" "${plas_scaffolds} scaffolds found via plasFlow"
			plasmidsFoundviaplasFlow=1
		else
			printf "%-20s: %-8s : %s\\n" "plasmid Assembly" "ALERT" "No plasmid scaffold found?"
			if [[ "${status}" == "SUCCESS" ]]; then
				status="WARNING"
			fi
		fi
	else
		printf "%-20s: %-8s : %s\\n" "plasmid Assembly" "SUCCESS" "No plasmid scaffold found using plasmidSpades"
	fi
elif [[ "${dec_family}" == "Enterobacteriaceae" ]]; then
	printf "%-20s: %-8s : %s\\n" "plasmid Assembly" "FAILED" "/plasFlow not found"
	status="FAILED"
fi

#Check short scaffolds reduction script
if [[ -s "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" ]]; then
	# Count the number of '>' still remaining after trimming the contig file
	full_longies=">"
	full_longies=$(grep -c ${full_longies} "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta")
	# Calculate the number of lost (short) scaffolds
	full_shorties=$(( full_scaffolds - full_longies ))
	if [ -z ${full_shorties} ]; then
		full_shorties=0
	fi
	#echo "${full_longies}"
	if [[ "${full_longies}" -le 200 ]]; then
		printf "%-20s: %-8s : %s\\n" "Contig Trim" "SUCCESS" "${full_longies} scaffolds remain. ${full_shorties} were removed due to shortness"
	else
		printf "%-20s: %-8s : %s\\n" "Contig Trim" "WARNING" "${full_longies} scaffolds remain which is high. ${full_shorties} were removed due to shortness"
		if [[ "${status}" == "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
			status="WARNING"
		fi
	fi
else
	printf "%-20s: %-8s : %s\\n" "Contig Trim" "FAILED" "/Assembly/${1}_scaffolds_trimmed.fasta not found"
	status="FAILED"
fi

#Check short scaffolds reduction script for plasmid assembly
#echo "${plasmidsFoundviaplasFlow}-Found?"
if [[ "${plasmidsFoundviaplasFlow}" -eq 1 ]]; then
	if [[ -s "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta" ]]; then
		# Count the number of '>' still remaining after trimming the contig file
		plas_longies=">"
		plas_longies=$(grep -c ${plas_longies} "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta")
		# Calculate the number of lost (short) scaffolds
		plas_shorties=$(( plas_scaffolds - plas_longies ))
		if [ -z ${plas_shorties} ]; then
			plas_shorties=0
		fi
		if [[ "${plas_longies}" -gt 0 ]]; then
			printf "%-20s: %-8s : %s\\n" "Plasmids contig Trim" "SUCCESS" "${plas_longies} scaffolds remain. ${plas_shorties} were removed due to shortness"
		else
			printf "%-20s: %-8s : %s\\n" "Plasmids contig Trim" "SUCCESS" "No plasmid scaffold found"
		fi
	elif [[ -f "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta" ]]; then
		printf "%-20s: %-8s : %s\\n" "Plasmids contig Trim" "SUCCESS" "No plasmid scaffolds found"
	else
		printf "%-20s: %-8s : %s\\n" "Plasmids contig Trim" "FAILED" "plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta not found"
		status="FAILED"
	fi
fi

#Check kraken on assembly
kraken_post_success=false
if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled.kraken" ]]; then
	#printf "%-20s: %-8s : %s\\n" "kraken postassembly" "SUCCESS" "Found"
	kraken_post_success=true
else
	printf "%-20s: %-8s : %s\\n" "kraken postassembly" "FAILED" "/kraken/postAssembly/${1}_paired.kraken not found"
	status="FAILED"
fi
#Check Krona output of assembly
if [[ "${kraken_post_success}" = true ]]; then
	if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled.krona" ]] && [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled.html" ]]; then
		#printf "%-20s: %-8s : %s\\n" "krona-kraken-pstasmb" "SUCCESS" "Found"
		:
	else
		printf "%-20s: %-8s : %s\\n" "krona-kraken-pstasmb" "FAILED" "/kraken/postAssembly/${1}_assembled.krona & /kraken/postAssembly/${1}_assembled.html not found"
		status="FAILED"
	fi
else
	printf "%-20s: %-8s : %s\\n" "krona-kraken-pstasmb" "FAILED" "preassembly kraken did not complete successfully"
	status="FAILED"
fi
#Check extraction and unclassified values for kraken post assembly
if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" ]]; then
	# Extracts many elements of the summary file to report unclassified and species classified reads and percentages
	unclass=$(head -n 1 "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" | cut -d' ' -f2)
	#true_unclass=$(head -n 1 "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
	domain=$(sed -n '2p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" | cut -d' ' -f2)
	genuspost=$(sed -n '7p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" | cut -d' ' -f4)
	speciespost=$(sed -n '8p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" | cut -d' ' -f4)
	speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" | cut -d' ' -f2)
	#true_speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
	# If there are no reads at the domain level, then report no classified reads
	if (( $(echo "${domain} <= 0" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "post Classify" "FAILED" "There are no classified reads (Did post assembly kraken fail too?)"
		status="FAILED"
	# If there are classified reads then check to see if percent unclassifed falls above the threshold limit. Report warning if too high or success and stats if below
	else
		if (( $(echo "${unclass} > ${unclass_flag}" | bc -l) )); then
			printf "%-20s: %-8s : %s\\n" "post Classify" "WARNING" "unclassified reads comprise ${unclass}% of total ${true_unclass}%"
			if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
				status="WARNING"
			fi
		elif (( $(echo "${speciespercent} < 50" | bc -l) )); then
			printf "%-20s: %-8s : %s\\n" "post Classify" "WARNING" "${genuspost} ${speciespost} is under 50% (${speciespercent}), possibly contaminated or contigs are weighted unevenly"
			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "SUCCESS" ]]; then
				status="WARNING"
			fi
		else
			#printf "%-20s: %-8s : %s\\n" "post Classify" "SUCCESS" "${speciespercent}%${true_speciespercent%} ${genuspost} ${speciespost} with ${unclass}%${true_unclass%} unclassified reads"
			printf "%-20s: %-8s : %s\\n" "post Classify" "SUCCESS" "${speciespercent}% ${genuspost} ${speciespost} with ${unclass}% unclassified reads"
		fi
	fi
# If no summary file was found
else
	printf "%-20s: %-8s : %s\\n" "post Classify" "FAILED" "/kraken/postAssembly/${1}_kraken_summary_assembled.txt not found"
	status="FAILED"
fi
#Check weighted kraken on assembly
kraken_weighted_success=false
if [[ ! -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled_BP.kraken" ]]; then
	if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled.kraken" ]]; then
		${shareScript}/run_kraken.sh "${1}" "post" "assembled" "${2}"
	fi
fi

# Quick separate check for contamination by finding # of species above ${contamination_threshold} in list file from kraken
if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled.list" ]]; then
	number_of_species=0
	while IFS= read -r line; do
		arrLine=(${line})
		# First element in array is the percent of reads identified as the current taxa
		percent=${arrLine[0]}
		percent_integer=$(echo "${percent}" | cut -d'.' -f1)
		# 3rd element is the taxon level classification
		classification=${arrLine[3]}
		#echo "${percent_integer} - ${contamination}"
		if [[ "${classification}" == "S" ]] && (( percent_integer > contamination_threshold )); then
			number_of_species=$(( number_of_species + 1 ))
		fi
	done < ${OUTDATADIR}/kraken/postAssembly/${1}_assembled.list
	if [[ $number_of_species -gt 1 ]]; then
		printf "%-20s: %-8s : %s\\n" "post Class Contam." "ALERT" "${number_of_species} species have been found above the ${contamination_threshold}% threshold"
		if [[ "${status}" == "SUCCESS" ]]; then
			status="ALERT"
		fi
	elif [[ "${number_of_species}" -eq 1 ]]; then
		:
	else
		printf "%-20s: %-8s : %s\\n" "post Class Contam." "ALERT" "No species have been found above ${contamination_threshold}% abundance"
		if [[ "${status}" == "SUCCESS" ]]; then
			status="ALERT"
		fi
	fi
	#echo "Number of species: ${number_of_species}"
fi



if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled_BP.kraken" ]]; then
	#printf "%-20s: %-8s : %s\\n" "kraken weighted" "SUCCESS" "Found"
	kraken_weighted_success=true
else
	printf "%-20s: %-8s : %s\\n" "kraken weighted" "FAILED" "Top match is under 50%, likely contaminated"
	status="FAILED"
fi
#Check Krona output of weighted assembly
if [[ "${kraken_weighted_success}" = true ]]; then
	if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled_weighted.krona" ]] && [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled_weighted_BP_krona.html" ]]; then
		#printf "%-20s: %-8s : %s\\n" "krona-kraken-weight" "SUCCESS" "Found"
		:
	else
		printf "%-20s: %-8s : %s\\n" "krona-kraken-weight" "FAILED" "/kraken/postAssembly/${1}_assembled_weighted.krona & /kraken/postAssembly/${1}_assembled_weighted_BP_krona.html not found"
		status="FAILED"
	fi
else
	printf "%-20s: %-8s : %s\\n" "krona-kraken-weight" "FAILED" "weighted conversion analysis of assembly kraken did not complete successfully"
	status="FAILED"
fi
#Check extraction and unclassified values for weighted kraken post assembly
if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" ]]; then
	# Extracts many elements of the summary file to report unclassified and species classified reads and percentages
	unclass=$(head -n 1 "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" | cut -d' ' -f2)
	#true_unclass=$(head -n 1 "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
	domain=$(sed -n '2p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" | cut -d' ' -f2)
	genusweighted=$(sed -n '7p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" | cut -d' ' -f4)
	speciesweighted=$(sed -n '8p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" | cut -d' ' -f4)
	speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" | cut -d' ' -f2)
	#true_speciespercent=$(sed -n '8p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" | cut -d' ' -f3 | sed -r 's/[)]+/%)/g')
	# If there are no reads at the domain level, then report no classified reads
	if (( $(echo "${domain} <= 0" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "weighted Classify" "FAILED" "There are no classified reads (Did post assembly kraken fail too?)"
		status="FAILED"
	# If there are classified reads then check to see if percent unclassifed falls above the threshold limit. Report warning if too high or success and stats if below
	else
		if (( $(echo "${unclass} > ${unclass_flag}" | bc -l) )); then
			printf "%-20s: %-8s : %s\\n" "weighted Classify" "WARNING" "unclassified reads comprise ${unclass}% of total ${true_unclass}%"
			if [ "${status}" = "SUCCESS" ] || [ "${status}" = "ALERT" ]; then
				status="WARNING"
			fi
		elif (( $(echo "${speciespercent} < 50" | bc -l) )); then
			printf "%-20s: %-8s : %s\\n" "weighted Classify" "FAILED" "${genusweighted} ${speciesweighted} is under 50% (${speciespercent}), likely contaminated"
			status="FAILED"
		else
			#printf "%-20s: %-8s : %s\\n" "weighted Classify" "SUCCESS" "${speciespercent}%${true_speciespercent%} ${genusweighted} ${speciesweighted} with ${unclass}%${true_unclass%} unclassified reads"
			printf "%-20s: %-8s : %s\\n" "weighted Classify" "SUCCESS" "${speciespercent}% ${genusweighted} ${speciesweighted} with ${unclass}% unclassified reads"
		fi
	fi
# If no summary file was found
else
	printf "%-20s: %-8s : %s\\n" "weighted Classify" "FAILED" "/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt not found"
	status="FAILED"
fi

# Quick separate check for contamination by finding # of species above ${contamination_threshold} in list file from kraken
if [[ -s "${OUTDATADIR}/kraken/postAssembly/${1}_assembled_BP.list" ]]; then
	number_of_species=0
	while IFS= read -r line; do
		arrLine=(${line})
		# First element in array is the percent of reads identified as the current taxa
		percent=${arrLine[0]}
		percent_integer=$(echo "${percent}" | cut -d'.' -f1)
		# 3rd element is the taxon level classification
		classification=${arrLine[3]}
		if [[ "${classification}" == "S" ]] && (( percent_integer > contamination_threshold )); then
			#echo "Adding ${line} because its S and greater than ${contamination_threshold}... ${percent_integer}"
			number_of_species=$(( number_of_species + 1 ))
		fi
	done < ${OUTDATADIR}/kraken/postAssembly/${1}_assembled_BP.list
	if [[ $number_of_species -gt 1 ]]; then
		printf "%-20s: %-8s : %s\\n" "weighted Contam." "FAILED" "${number_of_species} species have been found above the ${contamination_threshold}% threshold"
		status="FAILED"
	elif [[ "${number_of_species}" -eq 1 ]]; then
		:
	else
		printf "%-20s: %-8s : %s\\n" "weighted Contam." "FAILED" "No species have been found above the ${contamination_threshold}% threshold"
	fi
	#echo "Number of species: ${number_of_species}"
fi



#Check QUAST
if [[ -s "${OUTDATADIR}/Assembly_Stats/${1}_report.tsv" ]]; then
	# Extract the useful bits and report (to compare to Toms)
	contig_num=$(sed -n '14p' "${OUTDATADIR}/Assembly_Stats/${1}_report.tsv"| sed -r 's/[\t]+/ /g' | cut -d' ' -f3 )
	ass_length=$(sed -n '16p' "${OUTDATADIR}/Assembly_Stats/${1}_report.tsv" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
	N50=$(sed -n '18p' "${OUTDATADIR}/Assembly_Stats/${1}_report.tsv"  | sed -r 's/[\t]+/ /g'| cut -d' ' -f2)
	GC_con=$(sed -n '17p' "${OUTDATADIR}/Assembly_Stats/${1}_report.tsv" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
	printf "%-20s: %-8s : %s\\n" "QUAST" "SUCCESS" "#-${contig_num} length-${ass_length} n50-${N50} %GC-${GC_con}"
else
	printf "%-20s: %-8s : %s\\n" "QUAST" "FAILED" "/Assembly_Stats/report.tsv does not exist"
	status="FAILED"
fi

#Check QUAST on plasmid Assembly
if [[ "${plasmidsFoundviaplasFlow}" -eq 1 ]]; then
	if [[ -s "${OUTDATADIR}/Assembly_Stats_plasFlow/${1}_report.tsv" ]]; then
		# Extract the useful bits and report (to compare to Toms)
		contig_num_plas=$(sed -n '14p' "${OUTDATADIR}/Assembly_Stats_plasFlow/${1}_report.tsv"| sed -r 's/[\t]+/ /g' | cut -d' ' -f3 )
		ass_length_plas=$(sed -n '16p' "${OUTDATADIR}/Assembly_Stats_plasFlow/${1}_report.tsv" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
		N50_plas=$(sed -n '18p' "${OUTDATADIR}/Assembly_Stats_plasFlow/${1}_report.tsv"  | sed -r 's/[\t]+/ /g'| cut -d' ' -f2)
		GC_con_plas=$(sed -n '17p' "${OUTDATADIR}/Assembly_Stats_plasFlow/${1}_report.tsv" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
		printf "%-20s: %-8s : %s\\n" "QUAST_plasFlow" "SUCCESS" "#-${contig_num_plas} length-${ass_length_plas} n50-${N50_plas} %GC-${GC_con_plas}"
	else
		printf "%-20s: %-8s : %s\\n" "QUAST_plasFlow" "FAILED" "/Assembly_Stats_plasFlow/report.tsv does not exist"
		status="FAILED"
	fi
fi

# Get determinde taxonomy
if [[ ! -s "${OUTDATADIR}/${1}.tax" ]]; then
	"${shareScript}/determine_taxID.sh" "${1}" "${2}"
fi

source_call=$(head -n1 "${OUTDATADIR}/${1}.tax")
while IFS= read -r line; do
	# Grab first letter of line (indicating taxonomic level)
	first=${line:0:1}
	# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
	if [ "${first}" = "s" ]
	then
		dec_species=$(echo "${line}" | awk -F ' ' '{print $2}')
	elif [ "${first}" = "G" ]
	then
		dec_genus=$(echo "${line}" | awk -F ' ' '{print $2}')
	elif [ "${first}" = "F" ]
	then
		dec_family=$(echo "${line}" | awk -F ' ' '{print $2}')
	fi
done < "${OUTDATADIR}/${1}.tax"

if [[ "$dec_genus" != "Not_assigned" ]] && [[ "$dec_species" != "Not_assigned" ]]; then
	printf "%-20s: %-8s : %s\\n" "Taxa" "SUCCESS" "${dec_genus} ${dec_species}"
elif [[ "$dec_genus" != "Not_assigned" ]]; then
	printf "%-20s: %-8s : %s\\n" "Taxa" "FAILED" "None of the classifiers completed successfully"
elif [[ "$dec_species" != "Not_assigned" ]]; then
	printf "%-20s: %-8s : %s\\n" "Taxa" "WARNING" "No Species was able to be determined"
fi

# Check Assembly ratio
declare -A mmb_bugs
while IFS= read -r bug_lines || [ -n "$bug_lines" ]; do
	#bug_genus=$(echo "${bug_lines}" | cut -d'	' -f1)
	#bug_species=$(echo "${bug_lines}" | cut -d'	' -f2)
	bug_info=$(echo "${bug_lines}" | cut -d'	' -f4-)
	bug_size=$(echo "${bug_lines}" | cut -d'	' -f6)
	#bug_name="${bug_genus:0:1}.${bug_species}"
	bug_name=$(echo "${bug_lines}" | cut -d'	' -f3)
	#echo "Should be adding ${bug_size} for ${bug_name}"
	mmb_bugs["${bug_name}"]="${bug_size}"
done < ${local_DBs}/MMB_Bugs.txt
genus_initial="${dec_genus:0:1}"
ass_ID="${genus_initial}.${dec_species}"
#echo "${mmb_bugs[@]}"
#echo "${ass_ID}"
if [[ ! -z "${mmb_bugs[${ass_ID}]}" ]]; then
	#echo "Found Bug in DB: ${ass_ID}-${mmb_bugs[${ass_ID}]}"
	ass_ratio=$(awk -v p="${ass_length}" -v q="${mmb_bugs[${ass_ID}]}" 'BEGIN{printf("%.2f",p/q)}')
	if (( $(echo "$ass_ratio > 1.2" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Assembly ratio" "FAILED" "Too large - ${ass_ratio}x against ${ass_ID}"
		status="FAILED"
	elif (( $(echo "$ass_ratio < 0.8" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Assembly ratio" "FAILED" "Too small - ${ass_ratio}x against ${ass_ID}"
		status="FAILED"
	else
		printf "%-20s: %-8s : %s\\n" "Assembly ratio" "SUCCESS" "${ass_ratio}x against ${ass_ID}"
	fi
else
	printf "%-20s: %-8s : %s\\n" "Assembly_Ratio" "WARNING" "${ass_ID} does not exist in the DB"
	if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
		status="WARNING"
	fi
fi

# check coverage
if [[ -s "${OUTDATADIR}/preQCcounts/${1}_counts.txt" ]]; then
	line=$(head -n 1 "${OUTDATADIR}/preQCcounts/${1}_counts.txt")
	IFS='	' read -r -a qcs <<< "${line}"
	read_qc_info=${qcs[@]:1}
	# Extract q30 reads from qcCounts to calculate average coverage as q30_reads/assembly_length
	q30_reads=$(echo "${read_qc_info}" | awk -F ' ' '{print $2}')
	# Change later to AWK as this wont work on ASPEN, but consolidate won't likely be run on cluster
	if [[ ${ass_length} -gt 0 ]] && [[ ${q30_reads} -gt 0 ]]; then
		avg_coverage=$(bc <<<"scale=2 ; ${q30_reads} / ${ass_length}")
	else
		avg_coverage=0
	fi
	reads_low=40
	reads_high=90
	#echo "raw-${avg_coverage}"
	if (( $(echo "${avg_coverage} > ${reads_low}" | bc -l) )) && (( $(echo "${avg_coverage} < ${reads_high}" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Raw coverage" "SUCCESS" "${avg_coverage}x coverage based on raw reads"
	elif (( $(echo "${avg_coverage} > ${reads_high}" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Raw coverage" "ALERT" "${avg_coverage}x coverage based on raw reads"
		if [[ "${status}" == "SUCCESS" ]]; then
			status="ALERT"
		fi
	elif (( $(echo "${avg_coverage} < ${reads_low}" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Raw coverage" "FAILED" "${avg_coverage}x coverage based on raw reads"
		status="FAILED"
	fi
fi
if [[ -s "${OUTDATADIR}/preQCcounts/${1}_trimmed_counts.txt" ]]; then
	line=$(head -n 1 "${OUTDATADIR}/preQCcounts/${1}_trimmed_counts.txt")
	IFS='	' read -r -a qcs <<< "${line}"
	read_qc_info=${qcs[@]:1}
	# Extract q30 reads from qcCounts to calculate average coverage as q30_reads/assembly_length
	q30_reads=$(echo "${read_qc_info}" | awk -F ' ' '{print $2}')
	# Change later to AWK as this wont work on ASPEN, but consolidate won't likely be run on cluster
	if [[ ${ass_length} -gt 0 ]] && [[ ${q30_reads} -gt 0 ]]; then
		avg_coverage=$(bc <<<"scale=2 ; ${q30_reads} / ${ass_length}")
	else
		avg_coverage=0
	fi
	#echo "trimmed-${avg_coverage}"
	if (( $(echo "${avg_coverage} > ${reads_low}" | bc -l) )) && (( $(echo "${avg_coverage} < ${reads_high}" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Trimmed coverage" "SUCCESS" "${avg_coverage}x coverage based on trimmed reads"
	elif (( $(echo "${avg_coverage} > ${reads_high}" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Trimmed coverage" "ALERT" "${avg_coverage}x coverage based on trimmed reads"
		if [[ "${status}" == "SUCCESS" ]]; then
			status="ALERT"
		fi
	elif (( $(echo "${avg_coverage} < ${reads_low}" | bc -l) )); then
		printf "%-20s: %-8s : %s\\n" "Trimmed coverage" "ALERT" "${avg_coverage}x coverage based on trimmed reads"
		if [[ "${status}" == "SUCCESS" ]]; then
			status="ALERT"
		fi
	fi
fi

# Check prokka
if [[ -s "${OUTDATADIR}/prokka/${1}_PROKKA.gbf" ]]; then
	# Counts the number of genes present in the file using the 'CDS' identifier
	genes="CDS"
	genes=$(grep -c ${genes} "${OUTDATADIR}/prokka/${1}_PROKKA.gbf")
	printf "%-20s: %-8s : %s\\n" "prokka" "SUCCESS" "${genes} genes found"
elif [[ -s "${OUTDATADIR}/prokka/${1}_PROKKA.gbk" ]]; then
	# Counts the number of genes present in the file using the 'CDS' identifier
	genes="CDS"
	genes=$(grep -c ${genes} "${OUTDATADIR}/prokka/${1}_PROKKA.gbk")
	printf "%-20s: %-8s : %s\\n" "prokka" "SUCCESS" "${genes} genes found"
else
	printf "%-20s: %-8s : %s\\n" "prokka" "FAILED" "/prokka/${1}_PROKKA.gbf not found"
	status="FAILED"
fi

#Check BUSCO
if [[ -s "${OUTDATADIR}/BUSCO/short_summary_${1}.txt" ]]; then
	# Reads each line of the busco output file to extract the 3 that contain summary data to report
	while IFS= read -r line; do
		# If the line contains info for found buscos, total buscos, or database info grab it
		if [[ "${line}" == *"Complete BUSCOs (C)"* ]]
		then
			#echo "C-"${line}
			found_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
		elif [[ "${line}" == *"Total BUSCO groups searched"* ]];
		then
			#echo "T-"${line}
			total_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
		elif [[ "${line}" == *"The lineage dataset is:"* ]];
		then
			#echo "L-"${line}
			db=$(echo "${line}" | awk -F ' ' '{print $6}')
		fi
	done < "${OUTDATADIR}/BUSCO/short_summary_${1}.txt"
	percent_BUSCO_present=$(bc<<<"${found_buscos}*100/${total_buscos}")
	if [[ "${percent_BUSCO_present}" -gt 90 ]]; then
		printf "%-20s: %-8s : %s\\n" "BUSCO" "SUCCESS" "${percent_BUSCO_present}% (${found_buscos}/${total_buscos}) against ${db}"
	else
		printf "%-20s: %-8s : %s\\n" "BUSCO" "FAILED" "${percent_BUSCO_present}% (${found_buscos}/${total_buscos}) against ${db}"
		status="FAILED"
	fi
# If the busco summary file does not exist
else
	printf "%-20s: %-8s : %s\\n" "BUSCO" "FAILED" "/BUSCO/short_summary_${1}.txt not found"
	status="FAILED"
fi
#Check ANI
ani_found=false
# Goes through each file in the ANI folder to find if any match the proper format of the ANI output (this keeps old version and formats from sneaking in, pretty much uselss at this point)
if [[ -s "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_${dec_genus,}).txt" ]]; then
	mv "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_${dec_genus,}).txt" "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_${dec_genus^}).txt"
		 #"${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_${dec_genus^}).txt"
fi
if [[ -f "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_All).txt" ]]; then
	#echo "ALL"
	ani_info=$(head -n 1 "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_All).txt")
	ani_found=true
elif [[ -f "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_${dec_genus^}).txt" ]]; then
	ani_info=$(head -n 1 "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_${dec_genus^}).txt")
	ani_found=true
else
	for file in "${OUTDATADIR}/ANI/"*
	do
		if [[ "${file}" == *"best_ANI_hits_ordered(${1}_vs_"* ]]; then
			filename=${file}
			#echo "${OUTDATADIR}"
			#echo "${file}"
			#echo "${dec_genus^}"
			ani_info=$(head -n 1 "${file}")
			ani_found=true
			break
		fi
	done
fi
# Checks to see if the match boolean was toggled, if so it extracts the best match info and database from the string
if [[ "${ani_found}" = true ]]; then
	genusDB=$(echo "${filename##*/}" | cut -d'_' -f6 | cut -d')' -f1)
	percent_match="${ani_info:0:2}"
	#echo "${percent_match--}"
	if [[ "${percent_match}" = "0." ]]; then
		printf "%-20s: %-8s : %s\\n" "ANI" "FAILED" "No assembly file to work with"
		status="FAILED"
	else
		if [[ "${percent_match}" -ge 95 ]]; then
			printf "%-20s: %-8s : %s\\n" "ANI" "SUCCESS" "${ani_info} against ${genusDB}"
		else
			printf "%-20s: %-8s : %s\\n" "ANI" "FAILED" "${percent_match}% is too low, ${ani_info}"
			status="FAILED"
		fi
	fi
else
	if [[ "${dec_genus}" == "" ]]; then
		printf "%-20s: %-8s : %s\\n" "ANI" "FAILED" "Determine_TAXID did not discover a genus. Cant ANI on all samples, yet"
		status="FAILED"
	elif [[ ! -d "${OUTDATADIR}/ANI/" ]]; then
		printf "%-20s: %-8s : %s\\n" "ANI" "FAILED" "/ANI/ does not exist"
		status="FAILED"
	else
		printf "%-20s: %-8s : %s\\n" "ANI" "FAILED" "Attempted to compare to ${dec_genus} but /ANI/ does not have a best_ANI_hits file"
		status="FAILED"
	fi
fi

#Check c-SSTAR
if [[ -d "${OUTDATADIR}/c-sstar/" ]]; then
	if [[ ! -z "${3}" ]]; then
	 gapping="${3}"
	else
	 gapping="gapped"
	fi
	if [[ ! -z "${4}" ]]; then
		sim="${4}"
	else
		sim="98"
	fi
	csstar_file=$(find ${OUTDATADIR}/c-sstar/${1}.ResGANNCBI*.${gapping}_${sim}_sstar_summary.txt -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
	if [[ -z "${csstar_file}" ]]; then
		printf "%-20s: %-8s : %s\\n" "c-SSTAR" "FAILED" "/c-sstar/ does not have an sstar_summary file"
		status="FAILED"
	else
		header=$(head -n1 "${csstar_file}")
		ResGANNCBI_DB=$(echo "${csstar_file}" | rev | cut -d'.' -f3 | rev)
		#echo "${ResGANNCBI_DB} = ${ResGANNCBI_srst2_filename} ?"
		if [[ ${header} = *"No anti-microbial genes were found"* ]]; then
			if [[ "${ResGANNCBI_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
				printf "%-20s: %-8s : %s\\n" "c-SSTAR" "ALERT" "Completed, but NO KNOWN AMR genes were found in ${ResGANNCBI_DB} (DB up to date, as of ${today})"
			else
				printf "%-20s: %-8s : %s\\n" "c-SSTAR" "ALERT" "Completed, but NO KNOWN AMR genes were found in ${ResGANNCBI_DB} (DB NOT up to date! Most current DB: ${ResGANNCBI_srst2_filename})"
			fi
		elif [[ ${header} = "No Assembly found to run c-sstar with" ]]; then
			printf "%-20s: %-8s : %s\\n" "c-SSTAR" "FAILED" "No Assembly file to run through c-sstar"
			status="FAILED"
		else
			amr_genes_found=$(wc -l "${csstar_file}" | cut -d' ' -f1)
			# Prints out the counts of AR gene hits
			if [[ "${ResGANNCBI_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
				printf "%-20s: %-8s : %s\\n" "c-SSTAR" "SUCCESS" "${amr_genes_found} genes found in ${ResGANNCBI_DB} (DB up to date, as of ${today})"
			else
				printf "%-20s: %-8s : %s\\n" "c-SSTAR" "ALERT" "${amr_genes_found} genes found in ${ResGANNCBI_DB} (DB NOT up to date, Most current DB: ${ResGANNCBI_srst2_filename})"
			fi
		fi
	fi
else
	printf "%-20s: %-8s : %s\\n" "c-SSTAR" "FAILED" "/c-sstar/ does not exist"
	status="FAILED"
fi

# #Check c-SSTAR on plasmid Assembly
if [[ "${plasmidsFoundviaplasFlow}" -eq 1 ]]; then
	#Check c-SSTAR
	if [[ -d "${OUTDATADIR}/c-sstar_plasFlow/" ]]; then
		if [[ ! -z "${3}" ]]; then
			gapping="${3}"
		else
			gapping="gapped"
		fi
		if [[ ! -z "${4}" ]]; then
			sim="${4}"
		else
			sim="40"
		fi
		csstar_plasFlow_file=$(find ${OUTDATADIR}/c-sstar_plasFlow/${1}.ResGANNCBI*.${gapping}_${sim}_sstar_summary.txt -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
		if [[ -z "${csstar_plasFlow_file}" ]]; then
			printf "%-20s: %-8s : %s\\n" "c-SSTAR_plasFlow" "FAILED" "/c-sstar_plasFlow/ does not have an sstar_summary file"
			echo "Looking for ${OUTDATADIR}/c-sstar_plasFlow/${1}.ResGANNCBI.${gapping}_${sim}_sstar_summary.txt"
			status="FAILED"
		else
			header=$(head -n1 "${csstar_plasFlow_file}")
			ResGANNCBI_DB=$(echo "${csstar_plasFlow_file}" | rev | cut -d'.' -f3 | rev)
			if [[ ${header} = *"No anti-microbial genes were found"* ]]; then
				if [[ "${ResGANNCBI_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
					printf "%-20s: %-8s : %s\\n" "c-SSTAR_plasFlow" "ALERT" "Completed, but NO KNOWN AMR genes present from ${ResGANNCBI_DB} (DB up to date, as of ${today})"
					if [[ "${status}" == "SUCCESS" ]]; then
						status="ALERT"
					fi
				else
					printf "%-20s: %-8s : %s\\n" "c-SSTAR_plasFlow" "ALERT" "Completed, but NO KNOWN AMR genes present from ${ResGANNCBI_DB} (DB NOT up to date! Most current DB: ${ResGANNCBI_srst2_filename})"
					if [[ "${status}" == "SUCCESS" ]]; then
						status="ALERT"
					fi
				fi
			else
				amr_genes_found=$(wc -l "${csstar_plasFlow_file}" | cut -d' ' -f1)
				# Prints out the counts of AR gene hits
				if [[ "${ResGANNCBI_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
					printf "%-20s: %-8s : %s\\n" "c-SSTAR_plasFlow" "SUCCESS" "${amr_genes_found} genes found in ${ResGANNCBI_DB} (%ID defaults to 40) (DB up to date, as of ${today})"
				else
					printf "%-20s: %-8s : %s\\n" "c-SSTAR_plasFlow" "ALERT" "${amr_genes_found} genes found in ${ResGANNCBI_DB} (%ID defaults to 40) (DB NOT up to date! Most current DB: ${ResGANNCBI_srst2_filename})"
					if [[ "${status}" == "SUCCESS" ]]; then
						status="ALERT"
					fi
				fi
			fi
		fi
	else
		printf "%-20s: %-8s : %s\\n" "c-sstar_plasFlow" "FAILED" "/c-sstar_plasFlow/ does not exist - BOOYA"
		status="FAILED"
	# Signals that the current sample is completed with verification
	fi
fi

#Check c-SSTAR
if [[ -d "${OUTDATADIR}/GAMA/" ]]; then
	GAMA_file=$(find ${OUTDATADIR}/GAMA -maxdepth 1 -type f -name "${1}.ResGANNCBI*.GAMA"   -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
	if [[ -z "${GAMA_file}" ]]; then
		printf "%-20s: %-8s : %s\\n" "GAMA" "FAILED" "/GAMA/ does not have a .GAMA file"
		status="FAILED"
	else
		ResGANNCBI_DB=$(echo "${GAMA_file}" | rev | cut -d'.' -f2 | rev)
		#echo "${ResGANNCBI_DB} = ${ResGANNCBI_srst2_filename} ?"
		amr_genes_found=$(wc -l "${GAMA_file}" | cut -d' ' -f1)
		amr_genes_found=$(( amr_genes_found - 1))
		if [[ ${amr_genes_found} -le 0 ]]; then
			if [[ "${ResGANNCBI_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
				printf "%-20s: %-8s : %s\\n" "GAMA" "ALERT" "Completed, but NO KNOWN AMR genes were found in ${ResGANNCBI_DB} (DB up to date, as of ${today})"
			else
				printf "%-20s: %-8s : %s\\n" "GAMA" "ALERT" "Completed, but NO KNOWN AMR genes were found in ${ResGANNCBI_DB} (DB NOT up to date! Most current DB: ${ResGANNCBI_srst2_filename})"
			fi
		else
			# Prints out the counts of AR gene hits
			if [[ "${ResGANNCBI_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
				printf "%-20s: %-8s : %s\\n" "GAMA" "SUCCESS" "${amr_genes_found} genes found in ${ResGANNCBI_DB} (DB up to date, as of ${today})"
			else
				printf "%-20s: %-8s : %s\\n" "GAMA" "ALERT" "${amr_genes_found} genes found in ${ResGANNCBI_DB} (DB NOT up to date, Most current DB: ${ResGANNCBI_srst2_filename})"
			fi
		fi
	fi
else
	printf "%-20s: %-8s : %s\\n" "GAMA" "FAILED" "/GAMA/ does not exist"
	status="FAILED"
fi

# #Check c-SSTAR on plasmid Assembly
if [[ "${plasmidsFoundviaplasFlow}" -eq 1 ]]; then
	if [[ -d  "${OUTDATADIR}/GAMA_plasFlow" ]]; then
		#Check c-SSTAR
		GAMA_plasFlow_file=$(find ${OUTDATADIR}/GAMA_plasFlow -maxdepth 1 -type f -name "${1}.ResGANNCBI*.GAMA"   -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
		if [[ -z "${GAMA_plasFlow_file}" ]]; then
			printf "%-20s: %-8s : %s\\n" "GAMA_plasFlow" "FAILED" "/GAMA_plasFlow/ does not have a .GAMA file"
			status="FAILED"
		else
			ResGANNCBI_DB=$(echo "${GAMA_plasFlow_file}" | rev | cut -d'.' -f2 | rev)
			#echo "${ResGANNCBI_DB} = ${ResGANNCBI_srst2_filename} ?"
			plasmid_amr_genes_found=$(wc -l "${GAMA_plasFlow_file}" | cut -d' ' -f1)
			plasmid_amr_genes_found=$(( plasmid_amr_genes_found - 1))
			if [[ ${plasmid_amr_genes_found} -le 0 ]]; then
				if [[ "${ResGANNCBI_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
					printf "%-20s: %-8s : %s\\n" "GAMA_plasFlow" "ALERT" "Completed, but NO KNOWN AMR genes were found in ${ResGANNCBI_DB} (DB up to date, as of ${today})"
				else
					printf "%-20s: %-8s : %s\\n" "GAMA_plasFlow" "ALERT" "Completed, but NO KNOWN AMR genes were found in ${ResGANNCBI_DB} (DB NOT up to date! Most current DB: ${ResGANNCBI_srst2_filename})"
				fi
			else
				# Prints out the counts of AR gene hits
				if [[ "${ResGANNCBI_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
					printf "%-20s: %-8s : %s\\n" "GAMA_plasFlow" "SUCCESS" "${plasmid_amr_genes_found} genes found in ${ResGANNCBI_DB} (DB up to date, as of ${today})"
				else
					printf "%-20s: %-8s : %s\\n" "GAMA_plasFlow" "ALERT" "${plasmid_amr_genes_found} genes found in ${ResGANNCBI_DB} (DB NOT up to date, Most current DB: ${ResGANNCBI_srst2_filename})"
				fi
			fi
		fi
	else
		printf "%-20s: %-8s : %s\\n" "GAMA_plasFlow" "FAILED" "/GAMA_plasFlow/ does not exist"
		status="FAILED"
	fi
fi

# check SRST2 output
if [[ -d "${OUTDATADIR}/srst2/" ]]; then
	ResGANNCBI_srst2_file=$(find ${OUTDATADIR}/srst2/${1}__genes__ResGANNCBI*_srst2__results.txt -maxdepth 1 -type f -printf '%p\n' | sort -k6,6 -rt '_' -n | head -n 1)
	#echo ${ResGANNCBI_srst2_file}
	if [[ -s "${ResGANNCBI_srst2_file}" ]]; then
		ResGANNCBI_srst2_DB=$(echo "${ResGANNCBI_srst2_file}" | rev | cut -d'_' -f4,5 | rev)
		info_ResGANNCBI_List=$(head -n 1 "${ResGANNCBI_srst2_file}")
		IFS='	' read -r -a ResGANNCBI_array <<< "${info_ResGANNCBI_List}"
		ResGANNCBI_Num="${#ResGANNCBI_array[@]}"
		ResGANNCBI_Num=$(( ResGANNCBI_Num - 1 ))
		#echo "${info_ResGANNCBI_List} - ${ResGANNCBI_Num}"
		if [[ "${ResGANNCBI_Num}" -eq 0 ]]; then
			if [[ "${ResGANNCBI_srst2_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
				printf "%-20s: %-8s : %s\\n" "srst2" "ALERT" "Completed, but NO KNOWN AMR genes present from ${ResGANNCBI_srst2_DB} (DB up to date, as of ${today})"
				if [[ "${status}" == "SUCCESS" ]]; then
					status="ALERT"
				fi
			else
				printf "%-20s: %-8s : %s\\n" "srst2" "ALERT" "Completed, but NO KNOWN AMR genes present from ${ResGANNCBI_srst2_DB} (DB NOT up to date! Most current DB: ${ResGANNCBI_srst2_filename})"
				if [[ "${status}" == "SUCCESS" ]]; then
					status="ALERT"
				fi
			fi
		else
			if [[ "${ResGANNCBI_srst2_DB}" = "${ResGANNCBI_srst2_filename}" ]]; then
				printf "%-20s: %-8s : %s\\n" "srst2" "SUCCESS" "${ResGANNCBI_Num} genes found in ${ResGANNCBI_srst2_DB} (DB up to date, as of ${today})"
			else
				printf "%-20s: %-8s : %s\\n" "srst2" "ALERT" "${ResGANNCBI_Num} genes found in ${ResGANNCBI_srst2_DB} (DB NOT up to date! Most current DB: ${ResGANNCBI_srst2_filename})"
				if [[ "${status}" == "SUCCESS" ]]; then
					status="ALERT"
				fi
			fi
		fi
	else
		printf "%-20s: %-8s : %s\\n" "srst2" "FAILED" "genes file does not exist"
	fi
else
	printf "%-20s: %-8s : %s\\n" "srst2" "FAILED" "/srst2/ does not exist"
fi

# check MLST
if [[ -d "${OUTDATADIR}/MLST/" ]]; then
	if [[ -s "${OUTDATADIR}/MLST/${1}_Pasteur.mlst" ]] || [[ -s "${OUTDATADIR}/MLST/${1}.mlst" ]]; then
		if [[ -f "${OUTDATADIR}/MLST/${1}.mlst" ]]; then
			mv "${OUTDATADIR}/MLST/${1}.mlst" "${OUTDATADIR}/MLST/${1}_Pasteur.mlst"
		fi
		if [[ -f "${OUTDATADIR}/MLST/${1}_ecoli_2.mlst" ]]; then
			mv "${OUTDATADIR}/MLST/${1}_Pasteur.mlst" "${OUTDATADIR}/MLST/${1}_Achtman.mlst"
			mv "${OUTDATADIR}/MLST/${1}_ecoli_2.mlst" "${OUTDATADIR}/MLST/${1}_Pasteur.mlst"
		fi
		info=$(head -n 1 "${OUTDATADIR}/MLST/${1}_Pasteur.mlst")
		mlstype=$(echo "${info}" | cut -d'	' -f3)
		mlstdb=$(echo "${info}" | cut -d'	' -f2)
		#echo "'${mlstdb}:${mlstype}'"
		mlstdb="${mlstdb}(Pasteur)"
		if [ "${mlstdb}" = "-" ]; then
			if [ "${dec_genus}" ] && [ "${dec_species}" ]; then
				printf "%-20s: %-8s : %s\\n" "MLST" "WARNING" "no scheme found, check pubmlst for ${dec_genus} ${dec_species}"
				if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
					status="WARNING"
				fi
			else
				printf "%-20s: %-8s : %s\\n" "MLST" "FAILED" "no scheme found, check upstream as no genus has been assigned"
			fi
		elif [ "${mlstype}" = "-" ] || [ "${mlstype}" = "SUB" ]; then
			printf "%-20s: %-8s : %s\\n" "MLST" "WARNING" "no type found, possibly new type? Adding to maintenance_To_Do list"
			report_info=$(echo "${info}" | cut -d' ' -f2-)
			echo "${2}/${1}: Possible new MLST type - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
				status="WARNING"
			fi
		elif [ "${mlstype}" = "AU" ]; then
			printf "%-20s: %-8s : %s\\n" "MLST" "FAILED" "1+ allele is missing, cant determine ST type in ${mlstdb}"
			status="FAILED"
		else
			printf "%-20s: %-8s : %s\\n" "MLST" "SUCCESS" "TYPE is ${mlstype} from ${mlstdb}"
		fi
	else
		printf "%-20s: %-8s : %s\\n" "MLST" "FAILED" "${1}.mlst does not exist"
		status="FAILED"
	fi
	if [[ "${dec_genus}" = "Acinetobacter" ]]; then
		if [[ -s "${OUTDATADIR}/MLST/${1}_abaumannii.mlst" ]] || [[ -s "${OUTDATADIR}/MLST/${1}_Oxford.mlst" ]]; then
			if [[ -s "${OUTDATADIR}/MLST/${1}_abaumannii.mlst" ]]; then
				mv "${OUTDATADIR}/MLST/${1}_abaumannii.mlst" "${OUTDATADIR}/MLST/${1}_Oxford.mlst"
			fi
			info=$(tail -n 1 "${OUTDATADIR}/MLST/${1}_Oxford.mlst")
			mlstype=$(echo "${info}" | cut -d'	' -f3)
			mlstdb=$(echo "${info}" | cut -d'	' -f2)
			#echo "'${mlstdb}:${mlstype}'"
			if [ "${mlstdb}" = "abaumannii" ]; then
				if [ "${mlstype}" = "SUB" ] || [ "${mlstype}" = "-" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST" "WARNING" "no type found, possibly new type? Adding to maintenance_To_Do list"
					report_info=$(echo "${info}" | cut -d' ' -f2-)
					echo "${2}/${1}: Possible new MLST type - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
					if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
						status="WARNING"
					fi
				elif [ "${mlstype}" = "AU" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST" "FAILED" "1+ allele is missing, cant determine ST type in ${mlstdb}(Oxford)"
					status="FAILED"
				else
					printf "%-20s: %-8s : %s\\n" "MLST" "SUCCESS" "TYPE is ${mlstype} from ${mlstdb}(Oxford)"
				fi
			else
				echo "Not reporting as name and analyis expected do not match"
			fi
		fi
	fi
	if [[ "${dec_genus}" = "Escherichia" ]]; then
		if [[ -s "${OUTDATADIR}/MLST/${1}_Achtman.mlst" ]]; then
			info=$(tail -n 1 "${OUTDATADIR}/MLST/${1}_Achtman.mlst")
			mlstype=$(echo "${info}" | cut -d'	' -f3)
			mlstdb=$(echo "${info}" | cut -d'	' -f2)
			#echo "'${mlstdb}:${mlstype}'"
			if [ "${mlstdb}" = "ecoli" ]; then
				if [ "${mlstype}" = "SUB" ] || [ "${mlstype}" = "-" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST" "WARNING" "no type found, possibly new type? Adding to maintenance_To_Do list"
					report_info=$(echo "${info}" | cut -d' ' -f2-)
					echo "${2}/${1}: Possible new MLST type - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
					if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
						status="WARNING"
					fi
				elif [ "${mlstype}" = "AU" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST" "FAILED" "1+ allele is missing, cant determine ST type in ${mlstdb}(Achtman)"
					status="FAILED"
				else
					printf "%-20s: %-8s : %s\\n" "MLST" "SUCCESS" "TYPE is ${mlstype} from ${mlstdb}(Achtman)"
				fi
			else
				echo "Not reporting as name and analyis expected do not match"
			fi
		fi
	fi

	# Check srst2 MLSTs
	num_srst2_mlsts=$(find ${OUTDATADIR}/MLST -type f -name "*_srst2_*.mlst" | wc -l)
	#echo "${num_srst2_mlsts}"
	if [[ "${num_srst2_mlsts}" -eq 0 ]]; then
		#echo "No mlst srst2 was attempted on this isolate (${1})"
		:
	elif [[ "${num_srst2_mlsts}" -eq 1 ]]; then
		srst_mlst=$(find ${OUTDATADIR}/MLST -type f -name "*_srst2_*.mlst")
		if [[ "${srst_mlst}" == *"-Standard.mlst" ]]; then
			new_srst_mlst=${srst_mlst/Standard/Pasteur}
			mv ${srst_mlst} ${new_srst_mlst}
			srst_mlst=${new_srst_mlst}
		fi
		mlstype=$(tail -n1 ${srst_mlst} | cut -d'	' -f2)
		mlstdb=$(echo "${srst_mlst}" | rev | cut -d'-' -f1 | cut -d'.' -f2 | rev )
		if [ "${mlstype}" = "SUB" ] || [ "${mlstype}" = "-" ]; then
			printf "%-20s: %-8s : %s\\n" "MLST-srst2" "WARNING" "no type found, possibly new type? Adding to maintenance_To_Do list"
			report_info=$(echo "${info}" | cut -d' ' -f2-)
			echo "${2}/${1}: Possible new MLST type - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
			if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
				status="WARNING"
			fi
		elif [ "${mlstype}" = "AU" ]; then
			printf "%-20s: %-8s : %s\\n" "MLST-srst2" "FAILED" "1+ allele is missing, cant determine ST type in ${mlstdb}"
			status="FAILED"
		else
			printf "%-20s: %-8s : %s\\n" "MLST-srst2" "SUCCESS" "TYPE is ${mlstype} from ${mlstdb}"
		fi
	elif [[ "${num_srst2_mlsts}" -eq 2 ]]; then
		if [[ "${dec_genus}" = "Acinetobacter" ]]; then
			if [[ -f "${OUTDATADIR}/MLST/${1}_srst2_Acinetobacter_baumannii#1-Oxford.mlst" ]]; then
				srst_mlst="${OUTDATADIR}/MLST/${1}_srst2_Acinetobacter_baumannii#1-Oxford.mlst"
				mlstype=$(tail -n1 ${srst_mlst} | cut -d'	' -f2)
				mlstdb="abaumannii(Oxford)"
				if [ "${mlstype}" = "SUB" ] || [ "${mlstype}" = "-" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "WARNING" "no type found, possibly new type? Adding to maintenance_To_Do list"
					report_info=$(echo "${info}" | cut -d' ' -f2-)
					echo "${2}/${1}: Possible new MLST type - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
					if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
						status="WARNING"
					fi
				elif [ "${mlstype}" = "AU" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "FAILED" "1+ allele is missing, cant determine ST type in ${mlstdb}"
					status="FAILED"
				else
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "SUCCESS" "TYPE is ${mlstype} from ${mlstdb}"
				fi
			fi
			if [[ -f "${OUTDATADIR}/MLST/${1}_srst2_Acinetobacter_baumannii#2-Pasteur.mlst" ]]; then
				srst_mlst="${OUTDATADIR}/MLST/${1}_srst2_Acinetobacter_baumannii#2-Pasteur.mlst"
				mlstype=$(tail -n1 ${srst_mlst} | cut -d'	' -f2)
				mlstdb="abaumannii_2(Pasteur)"
				if [ "${mlstype}" = "SUB" ] || [ "${mlstype}" = "-" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "WARNING" "no type found, possibly new type? Adding to maintenance_To_Do list"
					report_info=$(echo "${info}" | cut -d' ' -f2-)
					echo "${2}/${1}: Possible new MLST type - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
					if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
						status="WARNING"
					fi
				elif [ "${mlstype}" = "AU" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "FAILED" "1+ allele is missing, cant determine ST type in ${mlstdb}"
					status="FAILED"
				else
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "SUCCESS" "TYPE is ${mlstype} from ${mlstdb}"
				fi
			fi
		elif [[ "${dec_genus}" = "Escherichia" ]]; then
			if [[ -f "${OUTDATADIR}/MLST/${1}_srst2_Escherichia_coli#1-Achtman.mlst" ]]; then
				srst_mlst="${OUTDATADIR}/MLST/${1}_srst2_Escherichia_coli#1-Achtman.mlst"
				mlstype=$(tail -n1 ${srst_mlst} | cut -d'	' -f2)
				mlstdb="ecoli(Achtman)"
				if [ "${mlstype}" = "SUB" ] || [ "${mlstype}" = "-" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "WARNING" "no type found, possibly new type? Adding to maintenance_To_Do list"
					report_info=$(echo "${info}" | cut -d' ' -f2-)
					echo "${2}/${1}: Possible new MLST type - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
					if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
						status="WARNING"
					fi
				elif [ "${mlstype}" = "AU" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "FAILED" "1+ allele is missing, cant determine ST type in ${mlstdb}"
					status="FAILED"
				else
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "SUCCESS" "TYPE is ${mlstype} from ${mlstdb}"
				fi
			fi
			if [[ -f "${OUTDATADIR}/MLST/${1}_srst2_Escherichia_coli#2-Pasteur.mlst" ]]; then
				srst_mlst="${OUTDATADIR}/MLST/${1}_srst2_Escherichia_coli#2-Pasteur.mlst"
				mlstype=$(tail -n1 ${srst_mlst} | cut -d'	' -f2)
				mlstdb="ecoli_2(Pasteur)"
				if [ "${mlstype}" = "SUB" ] || [ "${mlstype}" = "-" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "WARNING" "no type found, possibly new type? Adding to maintenance_To_Do list"
					report_info=$(echo "${info}" | cut -d' ' -f2-)
					echo "${2}/${1}: Possible new MLST type - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
					if [[ "${status}" = "SUCCESS" ]] || [[ "${status}" = "ALERT" ]]; then
						status="WARNING"
					fi
				elif [ "${mlstype}" = "AU" ]; then
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "FAILED" "1+ allele is missing, cant determine ST type in ${mlstdb}"
					status="FAILED"
				else
					printf "%-20s: %-8s : %s\\n" "MLST-srst2" "SUCCESS" "TYPE is ${mlstype} from ${mlstdb}"
				fi
			fi
		else
			printf "%-20s: %-8s : %s\\n" "MLST-srst2" "ALERT" "More than 1 srst file found for non AB or ecoli sample, look into it?"
			if [[ "${status}" == "SUCCESS" ]]; then
				status="ALERT"
			fi
		fi
	else
		printf "%-20s: %-8s : %s\\n" "MLST-srst2" "ALERT" "More than 2 srst files found, look into it?"
		if [[ "${status}" == "SUCCESS" ]]; then
			status="ALERT"
		fi
	fi
	# No MLST folder exists (pipeline must have failed as it would create a default one otherwise)
else
	printf "%-20s: %-8s : %s\\n" "MLST" "FAILED" "/MLST/ does not exist"
	status="FAILED"
fi
# check 16s Identification
if [[ -d "${OUTDATADIR}/16s/" ]]; then
	if [[ -s "${OUTDATADIR}/16s/${1}_16s_blast_id.txt" ]]; then
		info_b=$(head -n 1 "${OUTDATADIR}/16s/${1}_16s_blast_id.txt")
		genus_b=$(echo ${info_b} | cut -d' ' -f3)
		species_b=$(echo ${info_b} | cut -d' ' -f4-)
		IFS=' ' read -r -a id_array <<< "${info_b}"
		if [ ${#id_array[@]} -gt 3 ]; then
			extra_b="${id_array[@]:3:}"
		else
			extra_b=""
		fi
		#echo "g-${genus_b},s-${species_b}"
		if [ ! -z "${genus_b}" ] && [ ! -z "${species_b}" ]; then
			#if [[ "${genus_b}" == "No" ]] && [[ "${species_b}" == "16s" ]]; then
			#	printf "%-20s: %-8s : %s\\n" "16s_best_hit" "Warning" "No 16s sequences found"
			##		status="Warning"
			#	fi
			#else
				printf "%-20s: %-8s : %s\\n" "16s_best_hit" "SUCCESS" "${genus_b} ${species_b} ${extra_b}"
			#fi
		elif [ -z "${species_b}" ]; then
			if [ "${genus_b}" = "No_16s_sequences_found" ]; then
				printf "%-20s: %-8s : %s\\n" "16s_best_hit" "FAILED" "No 16s sequences found"
				status="FAILED"
			elif [ "${genus_b}" = "No_16s_matches_found" ]; then
				printf "%-20s: %-8s : %s\\n" "16s_best_hit" "FAILED" "16s sequences were found but were not able to be classified"
				status="FAILED"
			else
				printf "%-20s: %-8s : %s\\n" "16s_best_hit" "Warning" "Genus=${genus_b}, but no species found, Adding to maintenance_To_Do list"
				if [ "$status" = "SUCCESS" ]; then
					status="Warning"
				fi
			fi
			#report_info=$(echo "${info_b}" | cut -d' ' -f2-)
			#echo "${2}/${1}: 16s ID Warning - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"

		elif [ -z "${genus_b}" ]; then
			printf "%-20s: %-8s : %s\\n" "16s_best_hit" "FAILED" "No genus found, Adding to maintenance_To_Do list"
			report_info=$(echo "${info_b}" | cut -d' ' -f2-)
			echo "${2}/${1}: 16s ID Failure - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
			status="FAILED"
		else
			printf "%-20s: %-8s : %s\\n" "16s_best_hit" "FAILED" "Nothing found in ${1}_16s_blast_id.txt, Adding to maintenance_To_Do list"
			report_info=$(echo "${info_l}" | cut -d' ' -f2-)
			echo "${2}/${1}: 16s ID Failure - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
			status="FAILED"
		fi
		info_l=$(tail -n 1 "${OUTDATADIR}/16s/${1}_16s_blast_id.txt")
		genus_l=$(echo ${info_l} | cut -d' ' -f3)
		species_l=$(echo ${info_l} | cut -d' ' -f4-)
		IFS=' ' read -r -a id_array <<< "${info_l}"
		if [ ${#id_array[@]} -gt 3 ]; then
			extra_l="${id_array[@]:3:}"
		else
			extra_l=""
		fi
		if [ ! -z "${genus_l}" ] && [ ! -z "${species_l}" ]; then
			#if [[ "${genus_l}" == "No" ]] && [[ "${species_l}" == "16s" ]]; then
			#	printf "%-20s: %-8s : %s\\n" "16s_largest_hit" "Warning" "No 16s sequences found"
			#	if [ "$status" = "SUCCESS" ]; then
			#		status="Warning"
			#	fi
			#else
				printf "%-20s: %-8s : %s\\n" "16s_largest_hit" "SUCCESS" "${genus_l} ${species_l} ${extra_l}"
			#fi
		elif [ -z "${species_l}" ]; then
			if [ "${genus_l}" = "No_16s_sequences_found" ]; then
				printf "%-20s: %-8s : %s\\n" "16s_largest_hit" "FAILED" "No 16s sequences found"
				status="FAILED"
			elif [ "${genus_l}" = "No_16s_matches_found" ]; then
				printf "%-20s: %-8s : %s\\n" "16s_largest_hit" "FAILED" "16s sequences were found but were not able to be classified"
				status="FAILED"
			else
				printf "%-20s: %-8s : %s\\n" "16s_largest_hit" "Warning" "Genus=${genus_l}, but no species found, Adding to maintenance_To_Do list"
				if [ "$status" = "SUCCESS" ]; then
					status="Warning"
				fi
			#printf "%-20s: %-8s : %s\\n" "16s_largest_hit" "Warning" "Genus=${genus_l}, but no species found, Adding to maintenance_To_Do list"
			#report_info=$(echo "${info_l}" | cut -d' ' -f2-)
			#echo "${2}/${1}: 16s ID Warning - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
			fi
		elif [ -z "${genus_l}" ]; then
			printf "%-20s: %-8s : %s\\n" "16s_largest_hit" "FAILED" "no genus found, Adding to maintenance_To_Do list"
			report_info=$(echo "${info_l}" | cut -d' ' -f2-)
			echo "${2}/${1}: 16s ID Failure - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
			status="FAILED"
		else
			printf "%-20s: %-8s : %s\\n" "16s_largest_hit" "FAILED" "nothing found in ${1}_16s_blast_id.txt, Adding to maintenance_To_Do list"
			report_info=$(echo "${info_l}" | cut -d' ' -f2-)
			echo "${2}/${1}: 16s ID Failure - ${report_info}" >> "${shareScript}/maintenance_To_Do.txt"
			status="FAILED"
		fi
	else
		printf "%-20s: %-8s : %s\\n" "16s" "FAILED" "${1}_16s_blast_id.txt does not exist"
		status="FAILED"
	fi
# No 16s folder exists (pipeline must have failed as it would create a default one otherwise)
else
	printf "%-20s: %-8s : %s\\n" "16s" "FAILED" "/16s/ does not exist"
	status="FAILED"
fi

# check plasmids
if [[ -d "${OUTDATADIR}/plasmid/" ]]; then
	if [[ -s "${OUTDATADIR}/plasmid/${1}_results_table_summary.txt" ]]; then
		number_of_plasmids=0
		while read line_in; do
			line_in=$(echo ${line_in} | cut -d' ' -f1)
			if [[ "${line_in}" = "No" ]] || [[ "${line_in}" = "Enterococcus,Streptococcus,Staphylococcus" ]] || [[ "${line_in}" = "Enterobacteriaceae" ]] || [[ "${line_in}" = "Plasmid" ]]; then
				:
			else
				number_of_plasmids=$(( number_of_plasmids + 1 ))
			fi
		done < "${OUTDATADIR}/plasmid/${1}_results_table_summary.txt"
		printf "%-20s: %-8s : %s\\n" "plasmid" "SUCCESS" "${number_of_plasmids} replicons were found in the full scaffold"
	else
		printf "%-20s: %-8s : %s\\n" "plasmid" "FAILED" "results_table_summary.txt does not exist"
		status="FAILED"
	fi
# No plasmid folder exists
else
	printf "%-20s: %-8s : %s\\n" "plasmid" "FAILED" "/plasmid/ does not exist"
	status="FAILED"
fi

# # check plasmids (on plasmidAssembly)
if [[ "plasmidsFoundviaplasFlow" -eq 1 ]]; then
	if [[ -d "${OUTDATADIR}/plasmid_on_plasFlow/" ]]; then
		if [[ -s "${OUTDATADIR}/plasmid_on_plasFlow/${1}_results_table_summary.txt" ]]; then
			number_of_plasmids=0
			while read line_in; do
				line_in=$(echo ${line_in} | cut -d' ' -f1)
				if [[ "${line_in}" = "No" ]] || [[ "${line_in}" = "Enterococcus,Streptococcus,Staphylococcus" ]] || [[ "${line_in}" = "Enterobacteriaceae" ]] || [[ "${line_in}" = "Plasmid" ]]; then
					:
				else
					number_of_plasmids=$(( number_of_plasmids + 1 ))
				fi
			done < "${OUTDATADIR}/plasmid/${1}_results_table_summary.txt"
			printf "%-20s: %-8s : %s\\n" "plasmid-plasmidAsmb" "SUCCESS" "${number_of_plasmids} replicons were found in the plasmid scaffold"
		else
			printf "%-20s: %-8s : %s\\n" "plasmid-plasmidAsmb" "FAILED" "results_table_summary.txt does not exist"
			status="FAILED"
		fi
	# No plasmid folder exists
	else
		printf "%-20s: %-8s : %s\\n" "plasmid-plasmidAsmb" "FAILED" "/plasmid_on_plasFlow/ does not exist"
		status="FAILED"
	fi
fi
echo "---------- ${1} completed as ${status} ----------"

if [ "${status}" = "WARNING" ] || [ "${status}" = "FAILED" ]; then
	echo "${2}/${1}: ${status}" >> "${shareScript}/maintenance_To_Do.txt"
fi








#Script exited gracefully (unless something else inside failed)
exit 0
