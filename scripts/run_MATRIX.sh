#!/bin/bash -l

#
# Description: Pulls out MLST, AR genes, and plasmid repicons and creates a mashtree for the listed samples and consolidates them into one sheet when run from an alternate or old database
#
# Usage ./outbreak_analysis.sh -l path_to_list -t analysis_type (MATRIX|SNV|BOTH) -n analysis_identifier(e.g. outbreak identifier) [-g gapped|ungapped (analysis ran)] [-s identity 80|95|98|99|100] [-k clobberness keep|clobber] [-c path_to_config_file] [-d path_to_alt_DB] [-r alternate_REFSEQ_database] [-x path_to_crosswalk_file_for_known_controls] [-w do_crosswalk_but_with_main_list_file]
#
# Output location: Parameter

# v1.1.2 (08/26/2021)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "./run_MARIX.sh -r run_name [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?r:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		c)
			echo "Option -c triggered, argument = ${OPTARG}"
			config=${OPTARG};;
		r)
			echo "Option -l triggered, argument = ${OPTARG}"
			run_name=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi


if [[ -f "${config}" ]]; then
	echo "Loading special config file - ${config}"
	. "${config}"
else
	echo "Loading default config file"
	if [[ ! -f "./config.sh" ]]; then
		cp ./config_template.sh ./config.sh
	fi
	. ./config.sh
	cwd=$(pwd)
	config="${cwd}/config.sh"
fi

# Checks for proper argumentation
if [[ -z "${run_name}" ]]; then
	echo "run name is not set, exiting"
	exit 1
fi

# Creates the output directory if it does not exist
output_directory=${output_dir}/${run_name}
if [[ ! -d ${output_directory} ]]; then
	echo "Something went wrong, run folder does not exist"
	exit
else
	if [[ -f ${output_directory}/${run_name}_list.txt ]]; then
		list_file=${output_directory}/${run_name}_list.txt
	else
		echo "Something went wrong, there is no sample list file in the run directory"
		exit
	fi
	if [[ -f ${output_directory}/config_${run_name}.sh ]]; then
		config=${output_directory}/config_${run_name}.sh
		. ${config}
	fi
fi

# Setting default
sim=$csim
### project_parser_Percent_sim and project_parser_Percent_length are set in config file


# # Checks that value given for % Identity is one of the presets for csstar
# if [[ "${sim}" != 80 ]] && [[ "${sim}" != 95 ]] && [[ "${sim}" != 98 ]] && [[ "${sim}" != 99 ]] && [[ "${sim}" != 100 ]]; then
# 	echo "Identity is not one of the presets for csstar and therefore will fail, so using default 98"
# 	sim=98
# else
# 	# Overrule config level value for this one run
# 	project_parser_Percent_identity=${sim}
# fi

# Set database names to use
. "${src}/get_latest_DBs.sh" "${local_DBs}"
database_path=$(get_srst2)
database_and_version=$(get_srst2_filename)
AR_database_date=$(echo ${database_and_version} | rev | cut -d'_' -f1 | rev)
ani_database_path=$(get_ANI_REFSEQ)
ani_database_and_version=$(get_ANI_REFSEQ_Date)


echo "AR path and version - ${database_path} and ${database_and_version}"
echo "ANI path and version - ${ani_database_path} and ${ani_database_and_version}"

# Remove any pre-existing files from previous runs
if [[ -f ${output_directory}/${run_name}-mlst_summary.txt ]]; then
	rm ${output_directory}/${run_name}-mlst_summary.txt
fi
if [[ -f ${output_directory}/${run_name}-csstar_summary.txt ]]; then
	rm ${output_directory}/${run_name}-csstar_summary.txt
fi
if [[ -f ${output_directory}/${run_name}-plasmid_summary.txt ]]; then
	rm ${output_directory}/${run_name}-plasmid_summary.txt
fi
if [[ -f ${output_directory}/${run_name}_AR_plasmid_report.tsv ]]; then
	rm ${output_directory}/${run_name}_AR_plasmid_report.tsv
fi
if [[ -f ${output_directory}/${run_name}-sample_summary.txt ]]; then
	rm ${output_directory}/${run_name}-sample_summary.txt
fi
if [[ -f ${output_directory}/${run_name}-GAMMA_summary.txt ]]; then
	rm ${output_directory}/${run_name}-GAMMA_summary.txt
fi
if [[ -f ${output_directory}/${run_name}-GAMMA_rejects.txt ]]; then
	rm ${output_directory}/${run_name}-GAMMA_rejects.txt
fi
if [[ -f ${output_directory}/${run_name}-srst2.txt ]]; then
	rm ${output_directory}/${run_name}-srst2.txt
fi
if [[ -f ${output_directory}/${run_name}-srst2_rejects.txt ]]; then
	rm ${output_directory}/${run_name}-srst2_rejects.txt
fi

run_csstar="false"
run_srst2="false"
run_GAMMA="false"
run_ANI="false"
> "${output_directory}/${run_name}-csstar_todo.txt"
> "${output_directory}/${run_name}-srst2_todo.txt"
> "${output_directory}/${run_name}-GAMMA_todo.txt"
> "${output_directory}/${run_name}-ANI_todo.txt"
> "${output_directory}/${run_name}-csstar_rejects.txt"
> "${output_directory}/${run_name}-srst2_rejects.txt"
> "${output_directory}/${run_name}-GAMMA_rejects.txt"

# Check that each isolate has been compared to the newest ResGANNCBI DB file
while IFS= read -r line || [ -n "$line" ]; do
	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	SAMPLE_DATADIR="${output_directory}/${sample_name}"
	echo "checking for ${SAMPLE_DATADIR}/c-sstar/${sample_name}.${database_and_version}.${csstar_gapping}_${sim}_sstar_summary.txt"
	if [[ -s "${SAMPLE_DATADIR}/c-sstar/${sample_name}.${database_and_version}.${csstar_gapping}_${sim}_sstar_summary.txt" ]];
	then
		#echo "${project}/${sample_name} has newest ResGANNCBI for normal csstar already"
		:
	else
		echo "${run_name}/${sample_name} - ccstar needs to be run against ${database_and_version} at ${sim}"
		echo "${run_name}/${sample_name}" >> "${output_directory}/${run_name}-csstar_todo.txt"
		run_csstar="true"
	fi
	echo "checking for ${SAMPLE_DATADIR}/srst2/${sample_name}__genes__${database_and_version}_srst2__results.txt"
	if [[ -s ${SAMPLE_DATADIR}/FASTQs/${sample_name}_R1_001.fastq ]] && [[ -s ${SAMPLE_DATADIR}/FASTQs/${sample_name}_R1_001.fastq ]] || [[ -s ${SAMPLE_DATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz ]] && [[ -s ${SAMPLE_DATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
		#echo "FASTQs exist"
		if [[ -f "${SAMPLE_DATADIR}/srst2/${sample_name}__fullgenes__${database_and_version}_srst2__results.txt" ]] || [[ -f "${SAMPLE_DATADIR}/srst2/${sample_name}__genes__${database_and_version}_srst2__results.txt" ]]; then
			#echo "${project}/${sample_name} has newest ResGANNCBI for srst2 already"
			:
		else
			echo "${run_name}/${sample_name} - SRST2 needs to be run against ${database_and_version}"
			echo "${run_name}/${sample_name}" >> "${output_directory}/${run_name}-srst2_todo.txt"
			run_srst2="true"
		fi
	fi
	echo "checking for ${SAMPLE_DATADIR}/GAMMA/${sample_name}.${database_and_version}.gamma"
	if [[ -s "${SAMPLE_DATADIR}/GAMMA/${sample_name}.${database_and_version}.gamma" ]];
	then
		#echo "${project}/${sample_name} has newest ResGANNCBI for normal csstar already"
		:
	else
		echo "${run_name}/${sample_name} - GAMMA needs to be run against ${database_and_version}"
		echo "${run_name}/${sample_name}" >> "${output_directory}/${run_name}-GAMMA_todo.txt"
		run_GAMMA="true"
	fi
done < ${list_file}

while IFS= read -r line || [ -n "$line" ]; do
	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	SAMPLE_DATADIR="${output_directory}/${sample_name}"
	echo "checking for ${SAMPLE_DATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_REFSEQ_${ani_database_and_version}).txt"
	if [[ ! -s "${SAMPLE_DATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_REFSEQ_${ani_database_and_version}).txt" ]]; then
		echo "${run_name}/${sample_name} - ANI needs to be run against REFSEQ_${ani_database_and_version}"
		echo "${run_name}/${sample_name}" >> "${output_directory}/${run_name}-ANI_todo.txt"
		echo "${run_name}" >> "${output_directory}/${run_name}-sum_n_stats_todo.txt"
		run_ANI="true"
	fi
done < ${list_file}

### Needs update/change with modular code to be able to fill in missing gaps of tools
while IFS= read -r line || [ -n "$line" ]; do
	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	SAMPLE_DATADIR="${output_directory}/${sample_name}"
	# run csstar code
done < "${output_directory}/${run_name}-csstar_todo.txt"

while IFS= read -r line || [ -n "$line" ]; do
	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	SAMPLE_DATADIR="${output_directory}/${sample_name}"
	# run csstar code
done < "${output_directory}/${run_name}-GAMMA_todo.txt"

while IFS= read -r line || [ -n "$line" ]; do
	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	SAMPLE_DATADIR="${output_directory}/${sample_name}"
	# run csstar code
done < "${output_directory}/${run_name}-srst2_todo.txt"

while IFS= read -r line || [ -n "$line" ]; do
	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	SAMPLE_DATADIR="${output_directory}/${sample_name}"
	# run csstar code
done < "${output_directory}/${run_name}-ANI_todo.txt"

# Loop through and extracts and formats AR genes found in all isolates, as well as the primary MLST type and plasmid replicons. Each are output to separate files. Any AR genes that do not meet the length or % identity are copied to the rejects file.
while IFS= read -r line; do
	sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	SAMPLE_DATADIR="${output_directory}/${sample_name}"

	if [[ -f "${SAMPLE_DATADIR}/c-sstar/${sample_name}.${database_and_version}.${csstar_gapping}_${sim}_sstar_summary.txt" ]]; then
		csstar_list=""
		ARDB_full="${SAMPLE_DATADIR}/c-sstar/${sample_name}.${database_and_version}.${csstar_gapping}_${sim}_sstar_summary.txt"
		#echo "${ARDB_full}"
		# Extracts all AR genes from normal csstar output file and creates a lits of all genes that pass the filtering steps
		while IFS= read -r line; do
			# exit if no genes were found for the sample
			if [[ -z "${line}" ]] || [[ "${line}" == *"No anti-microbial genes were found"* ]]; then
				break
			fi
			IFS='	' read -r -a ar_line <<< "$line"
			percent_ID="${ar_line[6]}"
			percent_length="${ar_line[9]}"
			conferred=$(echo "${ar_line[1]}" | rev | cut -d'_' -f2- | rev)
			contig_number=$(echo "${ar_line[5]}" | rev | cut -d'_' -f3 | rev)
			gene="${ar_line[4]}"
			# Ensure that the gene passes % identity and % length threhsolds for reporting
			if [[ ${percent_length} -ge ${project_parser_Percent_length} ]] && [[ ${percent_ID} -ge ${project_parser_Percent_identity} ]] ; then
				if [[ -z "${csstar_list}" ]]; then
					#	echo "First csstar: ${gene}"
					csstar_list="${gene}(${conferred,,})[${percent_ID}NT/${percent_length}:#${contig_number}]"
				else
					if [[ ${csstar_list} == *"${gene}"* ]]; then
						#	echo "${gene} already found in ${csstar_list}"
						:
					else
						#	echo "${gene} not found in ${csstar_list}...adding it"
						csstar_list="${csstar_list},${gene}(${conferred,,})[${percent_ID}NT/${percent_length}:#${contig_number}]"
					fi
				fi
				# If length is less than predetermined minimum (90% right now) then the gene is added to a rejects list to show it was outside acceptable limits
			else
				echo -e "${project}\t${sample_name}\tfull_assembly\t${line}" >> ${output_directory}/${run_name}-csstar_rejects.txt
			fi
		done < ${ARDB_full}
		if [[ -z "${csstar_list}" ]]; then
			#echo "EMPTY-${project}	${sample_name}	No AR genes discovered"
			echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${run_name}-csstar_summary.txt
			csstar_list="No AR genes discovered"
		else
			#echo "OCCUPADO-${project}	${sample_name}	:${csstar_list}:"
			echo "${project}	${sample_name}	${csstar_list}" >> ${output_directory}/${run_name}-csstar_summary.txt
		fi
	else
		echo "${project}	${sample_name}	NO CURRENT FILE" >> ${output_directory}/${run_name}-csstar_summary.txt
		others=$(ls "${SAMPLE_DATADIR}/csstar")
		echo "${others}"
		echo "NOT FOUND - ${SAMPLE_DATADIR}/c-sstar/${sample_name}.${database_and_version}.${csstar_gapping}_${sim}_sstar_summary.txt"
		csstar_list="NO CURRENT FILE"
	fi

	#echo "^^^^^^^^^^^^^^^^^^^ ${SAMPLE_DATADIR}/GAMMA/${sample_name}.${database_and_version}.GAMMA"
	if [[ -f "${SAMPLE_DATADIR}/GAMMA/${sample_name}.${database_and_version}.gamma" ]]; then
    GAMMA_list=""
		GARDB_full="${SAMPLE_DATADIR}/GAMMA/${sample_name}.${database_and_version}.gamma"
		while IFS= read -r line; do
			# exit if no genes were found for the sample
			if [[ -z "${line}" ]]; then
				break
			elif [[ "${line}" == *"DB	Resistance	Gene_Family	Gene	Contig	Start"* ]]; then
				continue
			fi
			IFS='	' read -r -a ar_line <<< "$line"
			percent_BP_ID=$(echo "${ar_line[13]}" | awk '{ printf "%d", ($1*100) }' )
			percent_codon_ID=$(echo "${ar_line[12]}" | awk '{ printf "%d", ($1*100) }' )
			percent_length=$(echo "${ar_line[14]}" | awk '{ printf "%d", ($1*100) }' )
			conferred=$(echo "${ar_line[1]}" | rev | cut -d'_' -f2- | rev)
			contig_number=$(echo "${ar_line[4]}" | rev | cut -d'_' -f3 | rev)
			gene="${ar_line[3]}"
			# Ensure that the gene passes % identity and % length threhsolds for reporting
			if [[ ${percent_length} -ge ${project_parser_Percent_length} ]] && [[ ${percent_codon_ID} -ge ${project_parser_Percent_identity} ]]; then
				if [[ -z "${GAMMA_list}" ]]; then
				#	echo "First GAMMA: ${gene}"
					GAMMA_list="${gene,,}(${conferred,,})[${percent_BP_ID}NT/${percent_codon_ID}AA/${percent_length}:#${contig_number}]"
				else
					if [[ ${GAMMA_list} == *"${gene}"* ]]; then
					#	echo "${gene} already found in ${GAMMA_list}"
						:
					else
					#	echo "${gene} not found in ${GAMMA_list}...adding it"
						GAMMA_list="${GAMMA_list},${gene,,}(${conferred,,})[${percent_BP_ID}NT/${percent_codon_ID}AA/${percent_length}:#${contig_number}]"
					fi
				fi
			# If length is less than predetermined minimum (90% right now) then the gene is added to a rejects list to show it was outside acceptable limits
			else
				echo -e "${project}\t${sample_name}\tfull_assembly\t${line}" >> ${output_directory}/${run_name}-GAMMA_rejects.txt
			fi
		done < ${GARDB_full}
		if [[ -z "${GAMMA_list}" ]]; then
			echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${run_name}-GAMMA_summary.txt
			GAMMA_list="No AR genes discovered"
		else
			echo "${project}	${sample_name}	${GAMMA_list}" >> ${output_directory}/${run_name}-GAMMA_summary.txt
		fi
	else
		echo "${project}	${sample_name}	NO CURRENT FILE" >> ${output_directory}/${run_name}-GAMMA_summary.txt
		echo "NOT FOUND - ${SAMPLE_DATADIR}/GAMMA/${sample_name}.${database_and_version}.gamma"
		GAMMA_list="NO CURRENT FILE"
	fi

	# Adding in srst2 output in a similar fashion as to how the csstar genes are output to the file.
	if [[ -f "${SAMPLE_DATADIR}/srst2/${sample_name}__genes__${database_and_version}_srst2__results.txt" ]]; then
		srst2_results=""
		if [[ -s "${SAMPLE_DATADIR}/srst2/${sample_name}__fullgenes__${database_and_version}_srst2__results.txt" ]]; then
			while IFS= read -r line || [ -n "$line" ]; do
			#	echo "Start"
				gene=$(echo "${line}" | cut -d'	' -f3)
				#ODD WAY to do this right now, must look into later, but
				confers=$(echo "${line}" | cut -d'	' -f14 | cut -d';' -f3)
			#	echo "${gene}-${confers}"
				if [[ "${confers}" = "annotation" ]]; then
					continue
				fi
				if [[ -z "${confers}" ]]; then
					#if [[ -n ${gene} ]]; then
					#	if [[ "${gene,,}" == "agly_flqn" ]]; then
					#		confers="aminoglycoside_and_fluoroquinolone_resistance"
					#	elif [[ "${gene,,}" == "tetracenomycinc" ]]; then
					#		confers="tetracenomycinC_resistance"
					#	else
							confers=$(echo "${confers,,}" | rev | cut -d';' -f2 | rev)
					#	fi
					#fi
				fi
				confers=${confers,,//_resistance/}
				if [[ "${AR_database_date}" -ge 20210507 ]]; then
					allele=$(echo "${line}" | cut -d'	' -f4)
				else
					allele=$(echo "${line}" | cut -d'	' -f4  | rev | cut -d'_' -f2- | rev)
				fi
				if [[ "${allele}" = "Zn-dependent"* ]]; then
					allele="${allele}_hydrolase"
				fi
				coverage=$(echo "${line}" | cut -d'	' -f5)
				depth=$(echo "${line}" | cut -d'	' -f6)
				diffs=$(echo "${line}" | cut -d'	' -f7)
				if [[ ${diffs} == *"trunc"* ]]; then
					allele="TRUNC-${allele}"
				fi
				uncertainty=$(echo "${line}" | cut -d'	' -f8)
				divergence=$(echo "${line}" | cut -d'	' -f9)
				``
				length=$(echo "${line}" | cut -d'	' -f10)
				percent_length=$(echo "$coverage / 1" | bc)
				if [[ "${divergence}" = "0.0" ]]; then
					percent_ID=100
				else
					percent_ID=$(echo "100 - (($divergence + 1) / 1)" | bc)
				fi
			#	echo "${allele}/${coverage}/${depth}/${diffs}/${uncertainty}/${divergence}/${length}/${percent_ID}/${percent_length}"
			# Filter genes based on thresholds for length and percent identity
				if [[ "${percent_ID}" -ge ${project_parser_Percent_identity} ]] && [[ "${percent_length}" -ge ${project_parser_Percent_length} ]]; then
					info_line="${allele,,}(${confers,,})[${percent_ID,,}NT/${percent_length,,}]"
					 if [[ -z "${srst2_results}" ]]; then
					 	srst2_results=${info_line}
					 else
					 	srst2_results="${srst2_results},${info_line}"
					 fi
				else
					if [[ ${line} = "Sample	DB	gene"* ]]; then
						:
					else
						echo ${line} >> ${output_directory}/${run_name}-srst2_rejects.txt
					fi
				fi
			done < "${SAMPLE_DATADIR}/srst2/${sample_name}__fullgenes__${database_and_version}_srst2__results.txt"
			#echo "Test1"
			if [[ -z "${srst2_results}" ]]; then
				#echo "1"
				echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${run_name}-srst2.txt
				srst2_results="No AR genes discovered"
			else
				#echo "2"
				echo "${project}	${sample_name}	${srst2_results}" >> ${output_directory}/${run_name}-srst2.txt
			fi
		else
			echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${run_name}-srst2.txt
			srst2_results="No AR genes discovered"
		fi

	else
		#echo "3"
		echo "${project}	${sample_name}	NO CURRENT FILE" >> ${output_directory}/${run_name}-srst2.txt
		echo "NOT FOUND - ${SAMPLE_DATADIR}/srst2/${sample_name}__fullgenes__${database_and_version}_srst2__results.txt"
		srst2_results="NO CURRENT FILE"
	fi

	# Extracts taxonomic info
	if [[ ! -f "${SAMPLE_DATADIR}/${sample_name}.tax" ]]; then
		"${shareScript}/determine_taxID.sh" -n "${sample_name}" -p "${project}" -c "${config}"
	fi
	tax_file="${SAMPLE_DATADIR}/${sample_name}.tax"
	sed -i '/^$/d' "${SAMPLE_DATADIR}/${sample_name}.tax"
	tax_header=$(head -n1 "${SAMPLE_DATADIR}/${sample_name}.tax")
	echo "${tax_header}"
	taxonomy_source_type=$(echo "${tax_header}" | cut -d'(' -f2 | cut -d')' -f1)
	taxonomy_source=$(echo "${tax_header}" | cut -d'-' -f4-)
	#echo "Test-${tax_header};${taxonomy_source_type};${taxonomy_source}"

	#echo "Looking at ${SAMPLE_DATADIR}/${sample_name}.tax"
	genus=$(head -n7 "${SAMPLE_DATADIR}/${sample_name}.tax" | tail -n1 | cut -d'	' -f2)
	species=$(head -n8 "${SAMPLE_DATADIR}/${sample_name}.tax" | tail -n1 | cut -d'	' -f2)
	taxonomy="${genus} ${species}"
	if [[ "${taxonomy_source_type}" = "ANI_REFSEQ_UTD" ]]; then
		confidence_info=$(head -n1 "${taxonomy_source}")
	else
		taxonomy_source_type=$(echo "${taxonomy_source_type}" | cut -d'(' -f2 | cut -d')' -f1)
		confidence_percent=$(echo "${tax_header}" | cut -d'-' -f2)
		confidence_info="NO_ANI...${taxonomy_source_type}=${confidence_percent}"
	fi

	# Pulls MLST type for sample and adds it to the summary file
	if [[ -f "${SAMPLE_DATADIR}/MLST/${sample_name}_Pasteur.mlst" ]]; then
		mlst=$(head -n 1 ${SAMPLE_DATADIR}/MLST/${sample_name}_Pasteur.mlst)
		mlst=${mlst//,/\/}
		alleles=$(echo "${mlst}" | cut -d'	' -f4-)
		alleles=${alleles//,/\/}
		echo "${alleles}"
		alleles=${alleles//	/.}
		echo "${alleles}"
		#alleles=${alleles// /.}
		#echo "${alleles}"
		mlst=$(echo "${mlst}" | cut -d'	' -f3)
		if [[ "${mlst}" == *"SUB" ]] || [[ "${mlst}" == "AU" ]] || [[ "${mlst}" == "-" ]]; then
			# Since srst2 has #2 as Pasteurs we dont need to check taxonomy
			if [[ -s "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}-Pasteur.mlst" ]]; then
				mlst_srst2_line1=$(head -n 1 "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}-Pasteur.mlst")
				mlst_srst2_line2=$(tail -n 1 "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}-Pasteur.mlst")
			elif [[ -s "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#2-Pasteur.mlst" ]]; then
				mlst_srst2_line1=$(head -n 1 "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#2-Pasteur.mlst")
				mlst_srst2_line2=$(tail -n 1 "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#2-Pasteur.mlst")
			fi

			if [[ ! -z "${mlst_srst2_line1}" ]] && [[ ! -z "${mlst_srst2_line2}" ]]; then
				mlst2=$(echo "${mlst_srst2_line2}" | cut -d'	' -f2)
				if [[ "${mlst2}" = *"?"* ]] || [[ "${mlst2}" = *"*"* ]] || [[ "${mlst2}" = *"-"* ]]; then
					echo "srst2 MLST is not perfect, will not use (${mlst2})"
				else
					mlst="^ST${mlst2}"
					alleles=""
					if [[ ${taxonomy} = "Escherichia coli" ]]; then
						for index in {3..10}; do
							current_name=$(echo "${mlst_srst2_line1}" | cut -d'	' -f${index})
							current_allele=$(echo "${mlst_srst2_line2}" | cut -d'	' -f${index})
							if [[ ${index} -eq 3 ]]; then
								alleles="${current_name}(${current_allele})"
							else
								alleles="${alleles}.${current_name}(${current_allele})"
							fi
						done
					else
						for index in {3..9}; do
							current_name=$(echo "${mlst_srst2_line1}" | cut -d'	' -f${index})
							current_allele=$(echo "${mlst_srst2_line2}" | cut -d'	' -f${index})
							if [[ ${index} -eq 3 ]]; then
								alleles="${current_name}(${current_allele})"
							else
								alleles="${alleles}.${current_name}(${current_allele})"
							fi
						done
					fi
				fi
			fi
		else
			mlst="ST${mlst}"
		fi
	else
		mlst="N/A"
		alleles="N/A"
	fi
	echo -e "${project}\t${sample_name}\t${mlst}\t${alleles}" >> ${output_directory}/${run_name}-mlst_summary.txt

	# Pulls Alternate MLST type for sample and adds it to the summary file
	if [[ -f "${SAMPLE_DATADIR}/MLST/${sample_name}_Oxford.mlst" ]]; then
		alt_mlst_file="${SAMPLE_DATADIR}/MLST/${sample_name}_Oxford.mlst"
	elif [[ -f "${SAMPLE_DATADIR}/MLST/${sample_name}_Achtman.mlst" ]]; then
		alt_mlst_file="${SAMPLE_DATADIR}/MLST/${sample_name}_Achtman.mlst"
	else
		alt_mlst_file=""
	fi
	if [[ -n "${alt_mlst_file}" ]]; then
		alt_mlst=$(tail -n 1 "${alt_mlst_file}")
		alt_alleles=$(echo "${alt_mlst}" | cut -d'	' -f4-)
		alt_alleles=${alt_alleles//	/.}
		alt_mlst=$(echo "${alt_mlst}" | cut -d'	' -f3)
		alt_mlst=${alt_mlst//,/\/}
		alt_alleles=${alt_alleles//,/\/}
		if [[ "${alt_mlst}" == *"SUB" ]] || [[ "${alt_mlst}" == "AU" ]] || [[ "${alt_mlst}" == "-" ]]; then
			# Since srst2 has #2 as Pasteurs we dont need to check taxonomy
			if [[ -s "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#1-Oxford.mlst" ]]; then
				mlst_srst2_line1=$(head -n 1 "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#1-Oxford.mlst")
				mlst_srst2_line2=$(tail -n 1 "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#1_Oxford.mlst")
			elif [[ -s "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#1-Achtman.mlst" ]]; then
				mlst_srst2_line1=$(head -n 1 "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#1-Achtman.mlst")
				mlst_srst2_line2=$(tail -n 1 "${SAMPLE_DATADIR}/MLST/${sample_name}_srst2_${taxonomy// /_}#1-Achtman.mlst")
			fi

			if [[ ! -z "${mlst_srst2_line1}" ]] && [[ ! -z "${mlst_srst2_line2}" ]]; then
				alt_mlst2=$(echo "${mlst_srst2_line2}" | cut -d'	' -f2)
				if [[ "${alt_mlst2}" = *"?"* ]] || [[ "${alt_mlst2}" = *"*"* ]] || [[ "${alt_mlst2}" = *"-"* ]]; then
					echo "srst2 MLST is not perfect, will not use (${alt_mlst2})"
				else
					alt_mlst="^ST${alt_mlst2}"
					alleles=""
					if [[ ${taxonomy} = "Escherichia coli" ]]; then
						for index in {3..10}; do
							current_name=$(echo "${mlst_srst2_line1}" | cut -d'	' -f${index})
							current_allele=$(echo "${mlst_srst2_line2}" | cut -d'	' -f${index})
							if [[ ${index} -eq 3 ]]; then
								alleles="${current_name}(${current_allele})"
							else
								alleles="${alleles}.${current_name}(${current_allele})"
							fi
						done
					else
						for index in {3..9}; do
							current_name=$(echo "${mlst_srst2_line1}" | cut -d'	' -f${index})
							current_allele=$(echo "${mlst_srst2_line2}" | cut -d'	' -f${index})
							if [[ ${index} -eq 3 ]]; then
								alleles="${current_name}(${current_allele})"
							else
								alleles="${alleles}.${current_name}(${current_allele})"
							fi
						done
					fi
				fi
			fi
		else
			alt_mlst="ST${alt_mlst}"
		fi
	else
		alt_mlst="N/A"
		alt_alleles="N/A"
	fi
	echo -e "${project}\t${sample_name}\t${alt_mlst}\t${alt_alleles}" >> ${output_directory}/${run_name}-alt_mlst_summary.txt

	# Print all extracted info to primary file
	echo -e "${project}\t${sample_name}\t${output_dir}\t${taxonomy}\t${taxonomy_source_type}\t${confidence_info}\t${mlst}\t${alleles}\t${alt_mlst}\t${alt_alleles}\t${csstar_list}\t${srst2_results}\t${GAMMA_list}" >> ${output_directory}/${run_name}-sample_summary.txt

	# Goes through the plasmid file of the sample and adds all found plasmid replicons to the summary file
	#echo "Starting plasmid extraction"
	if [[ -f ${SAMPLE_DATADIR}/plasmidFinder/${sample_name}_results_table_summary.txt ]]; then
		#echo "Found plasmid file"
		:
	fi
	full_contigs=">"
	full_contigs=$(grep -c ${full_contigs} "${SAMPLE_DATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta")
	added=0
	while IFS= read -r plasmid || [ -n "$plasmid" ]; do
		line_in=$(echo ${plasmid} | cut -d' ' -f1)
		if [[ "${line_in}" = "No" ]] || [[ "${line_in}" = "Enterococcus,Streptococcus,Staphylococcus" ]] || [[ "${line_in}" = "Enterobacteriaceae" ]] || [[ "${line_in}" = "Plasmid" ]]; then
			# echo "Not using line: $plasmid"
			:
		else
			echo -e "${project}\t${sample_name}\tfull_assembly\t${plasmid}" >> ${output_directory}/${run_name}-plasmid_summary.txt
			added=1
		fi
	done < ${SAMPLE_DATADIR}/plasmidFinder/${sample_name}_results_table_summary.txt
	if [[ "${added}" -eq 0 ]]; then
		echo -e "${project}\t${sample_name}\tfull_assembly\tNo_Plasmids_Found\t${full_contigs}_contigs-${components}_components" >> ${output_directory}/${run_name}-plasmid_summary.txt
	fi
done < ${list_file}

# Calls script that sorts and formats all isolates info into a matrix for easy viewing
python3 "${src}/matrix_maker.py" -s "${output_directory}/${run_name}-sample_summary.txt" -p "${output_directory}/${run_name}-plasmid_summary.txt" -o "${output_directory}/${run_name}_matrix.csv" -d "${database_and_version}" -m "${sim}" -l "${project_parser_Percent_length}"

if [[ ! -d "${output_directory}/matrix_files" ]]; then
	mkdir "${output_directory}/matrix_files"
fi
declare -a move_list
move_list=(csstar_todo GAMMA_todo srst2_todo alt_mlst_summary csstar_rejects csstar_summary GAMMA_rejects GAMMA_summary mlst_summary plasmid_summary sample_summary srst2 srst2_rejects ANI_todo)
for mlist in "${move_list[@]}"; do
	if [[ -f "${output_directory}/${run_name}-${mlist}.txt" ]]; then
 		mv "${output_directory}/${run_name}-${mlist}.txt" "${output_directory}/matrix_files/${run_name}-${mlist}.txt"
	else
		echo "${output_directory}/${run_name}-${mlist}.txt does not exist to move"
	fi
done

global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "run_MATRIX.sh for ${run_name} has completed" "${global_end_time}"
