#!/bin/bash -l

#$ -o getTax.out
#$ -e getTax.err
#$ -N getTax
#$ -cwd
#$ -q short.q

#
# Description: Creates a single file that attempts to pull the best taxonomic information from the isolate. Currently, it operates in a linear fashion, e.g. 1.ANI, 2.16s, 3.kraken, 4.Gottcha
# 	The taxon is chosen based on the highest ranked classifier first
#
# Usage: ./determine_texID.sh path_to_sample_folder path_to_databases
#
# Output: path_to_sample_folder/sample_name.tax
#
# Modules required: None
#
# v1.0.6b (05/12/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

### -------------------------------------------------------------------------- ###
#           REALLY NEED TO MAKE AN ARG PARSER IN HERE VERY SOON                  #
### -------------------------------------------------------------------------- ###


# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./determine_taxID.sh path_to_sample_folder path_to_databases"
	echo "Output is saved to path_to_sample_folder/sample_name.tax"
	exit 0
elif [[ ! -d "${1}" ]]; then
	echo "Path (${1}) does not exist, exiting determine_taxID.sh"
	exit 2
elif [[ -z "${2}" ]]; then
	echo "Empty sample_id supplied to determine_taxID.sh, exiting"
	exit 3
elif [[ ! -d "${2}" ]]; then
	echo "Database path (${2}) does not exist, exiting determine_taxID.sh"
	exit 4
elif [[ ! -f "${2}/taxes.csv" ]]; then
	echo "Database path (${2}) does not contain taxes.csv, upper classifcation can not be performed"
fi

OUTDATADIR="${1}"
databases="${2}"
# Based upon standard naming protocols pulling last portion of path off should result in proper sample name
sample_name=$(echo "${OUTDATADIR}" | rev | cut -d'/' -f1 | rev)
# Based upon standard naming protocols pulling 2nd to last portion of path off should result in proper project name
project_name=$(echo "${OUTDATADIR}" | rev | cut -d'/' -f2 | rev)


# Set default values for a ll taxonomic levels
Domain="Not_assigned"
Phylum="Not_assigned"
Class="Not_assigned"
Order="Not_assigned"
Family="Not_assigned"
Genus="Not_assigned"
species="Not_assigned"
source="Not_assigned"
confidence_index="0"
source_file="Not_assigned"


# Function to check which source to use as the 'determinator'. Single int parameter can be used to tell which level to jump in at
Check_source() {
	start_at="${1}"
	if [[ "${start_at}" -le 1 ]]; then
		for f in ${OUTDATADIR}/ANI/*; do
			#echo $f
			if [[ "${f}" = *"best_ANI_hits_ordered"* ]]; then
				header=$(head -n1 ${f})
				if [[ ${header} != "No matching ANI database found for"* ]] && [[ ${header} != "0.00%"* ]] ; then
		    	do_ANI
		    	return
				else
					Check_source 2
				fi
			fi
		done
	fi
	if [[ "${start_at}" -le 2 ]]; then
		if [[ -s "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt" ]]; then
			best_line=$(head -n1 "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt")
			largest_line=$(tail -n1 "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt")
			IFS='	' read -r -a best_array <<< "$best_line"
			IFS='	' read -r -a largest_array <<< "$largest_line"
			best_arr_size="${#best_array[@]}"
			largest_arr_size="${#largest_array[@]}"
			best_species=$(echo ${best_line} | cut -d'	' -f3)
			largest_species=$(echo ${largest_line} | cut -d'	' -f3)
			#echo "largest:${largest_species}:"
			#echo "best:${best_species}:"
			if [[ "${largest_arr_size}" -ge 3 ]]; then
				if [[ "${largest_array[2]}" == "Unidentified" ]] || [[ "${largest_array[2]}" == "No_16s_"* ]] || [[ "${largest_array[2]}" == "uncultured"* ]]; then
					:
				else
					do_16s "largest"
					return
				fi
			elif [[ "${best_arr_size}" -ge 3 ]] ; then
				if [[ "${best_array[2]}" == "Unidentified" ]]  || [[ "${best_array[2]}" == "No_16s_"* ]] || [[ "${best_array[2]}" == "uncultured"* ]]; then
					:
				else
					do_16s "best"
					return
				fi
			fi
		fi
	fi
	# if [[ "${start_at}" -le 3 ]];then
	# 	if [[ -s "${OUTDATADIR}/gottcha/${sample_name}_gottcha_species_summary.txt" ]]; then
	# 		do_GOTTCHA
	# 		return
	# 	fi
	# fi
	if [[ "${start_at}" -le 4 ]]; then
		if [[ -s "${OUTDATADIR}/kraken/postAssembly/${sample_name}_kraken_summary_assembled_BP.txt" ]]; then
			do_Kraken
		return
		fi
	fi
	echo "No ACCEPTABLE source found to determine taxonomy"
}

# Function to pull info from ANI output
do_ANI() {
	source="ANI"
	#echo "${source}"

  SCRIPT_DIR=$(dirname "$(readlink -f "$0")") # get dir of current script
  source ${SCRIPT_DIR}/get_latest_DBs.sh ${databases}
	refseq_ANI_date=$(get_ANI_REFSEQ_Date)
	#echo "${refseq_ANI_date}"
	if [[ -f "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_REFSEQ_${refseq_ANI_date}).txt" ]]; then
		source_file="${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_REFSEQ_${refseq_ANI_date}).txt"
		sample_name_underscore_count=$(tr -dc '_' <<<"${sample_name}" | awk '{ print length; }')
		reindex=$(( sample_name_underscore_count + 6))
		source="ANI_REFSEQ_UTD"
	else
		next_best_refseq_file=$(find ${OUTDATADIR}/ANI/ -maxdepth 1 -name "${sample_name}_vs_REFSEQ_*" -type f -printf '%p\n' | sort -k$reindex,$reindex -rt '_' -n | head -n 1)
		source_file=${next_best_refseq_file}
		source="ANI_REFSEQ_OOD"
	fi
	if [[ -f "${source_file}" ]]; then
		header=$(head -n 1 "${source_file}")
		#echo "${header}"
		Genus=$(echo "${header}" | cut -d' ' -f1 | cut -d'-' -f3)
		species=$(echo "${header}" | cut -d' ' -f2 | cut -d'(' -f1 | sed 's/[][]//g')
		confidence_index=$(echo "${header}" | cut -d' ' -f1 | cut -d'-' -f1,2)
		#echo "${Genus}-${species}"
	else
		Check_source 2
	fi
}

# Function to pull best info from 16s output (largest vs highest bit score)
do_16s() {
	source_file="${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt"
	if [[ "${1}" = "largest" ]]; then
		line=$(tail -n 1 "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt")
		source="16s_largest"
		if [[ -f "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all.sorted" ]]; then
			confidence_index=$(head -n1 "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all.sorted" | cut -d'	' -f3)
			confidence_index="${confidence_index}"
		else
			confidence_index=0
		fi
	elif [[ "${1}" = "best" ]]; then
		line=$(head -n 1 "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt")
		source="16s_best"
		if [[ -f "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all" ]]; then
			confidence_index=$(head -n1 "${OUTDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all" | cut -d'	' -f3)
			confidence_index="${confidence_index}"
		else
			confidence_index=0
		fi
	else
		break
	fi
	Genus=$(echo "${line}" | cut -d"	" -f3 | cut -d" " -f1)
	species=$(echo "${line}" | cut -d"	" -f3 | cut -d" " -f2)
}

# Function to pull info from gottcha output
do_GOTTCHA() {
	source="GOTTCHA"
	source_file="${OUTDATADIR}/gottcha/${sample_name}_gottcha_species_summary.txt"
	#echo "${source}"
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{print $5}')
		elif [ "${first}" = "G" ]
		then
			Genus=$(echo "${line}" | awk -F ' ' '{print $4}')
		fi
		confidence_index=$(tail -n1 "${source_file}" | cut -d' ' -f2)
		confidence_index="${confidence_index}"
	done < "${OUTDATADIR}/gottcha/${sample_name}_gottcha_species_summary.txt"
}

# Function to pull info from kraken output based on assembly
do_Kraken() {
	source="Kraken"
	source_file="${OUTDATADIR}/kraken/postAssembly/${sample_name}_kraken_summary_assembled_BP.txt"
	#echo "${source}"
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			species=$(echo "${line}" | awk -F ' ' '{print $4}')
		elif [ "${first}" = "G" ]
		then
			Genus=$(echo "${line}" | awk -F ' ' '{print $4}')
		fi
	done < "${OUTDATADIR}/kraken/postAssembly/${sample_name}_kraken_summary_assembled_BP.txt"
	confidence_index=$(tail -n1 "${source_file}" | cut -d' ' -f2)
	confidence_index="${confidence_index}"
}

# Start the program by checking ALL sources
Check_source 0

# Check if genus was assigned
if [[ ! -z ${Genus} ]]; then
	Genus=$(echo ${Genus} | tr -d [:space:] | tr -d "[]")
fi
# Check if species was assigned
if [[ ! -z ${species} ]]; then
	species=$(echo ${species} | tr -d [:space:])
fi

# Check if genus was assigned as peptoclostridium and relabel it as Clostridium for downstream analyses relying on this older naming convention
if [[ ${Genus} == "Peptoclostridium" ]]; then
	Genus="Clostridium"
fi

# Using premade database fill in upper levels of taxonomy info based on genus
while IFS= read -r line  || [ -n "$line" ]; do
	DB_genus=$(echo ${line} | cut -d"," -f1)
	#echo ":${Genus}:${DB_genus}:"
	if [[ "${Genus,}" = "${DB_genus}" ]]; then
			tax_DB="${databases}/taxes.csv"
			Domain=$(echo "${line}" | cut -d"," -f2)
			Phylum=$(echo "${line}" | cut -d"," -f3)
			Class=$(echo "${line}" | cut -d"," -f4)
			Order=$(echo "${line}" | cut -d"," -f5)
			Family=$(echo "${line}" | cut -d"," -f6 | tr -d '\r' )
			#echo ":${Family}:"
			break
	fi
done < "${databases}/taxes.csv"

# Print output to tax file for sample
echo -e "(${source})-${confidence_index}%-${source_file}\nD:	${Domain}\nP:	${Phylum}\nC:	${Class}\nO:	${Order}\nF:	${Family}\nG:	${Genus}\ns:	${species}\n" > "${OUTDATADIR}/${sample_name}.tax"
