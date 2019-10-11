#!/bin/sh -l

#$ -o kraken_weigh_contigs.out
#$ -e kraken_weigh_contigs.err
#$ -N kraken_weigh_contigs
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Grabs the best species match based on %/read hits from the kraken tool run
#
# Usage: ./kraken_weigh_contigs.sh sample_name run_ID kraken_version[kraken|kraken2]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/kraken(2)/
#
# No modules required
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [ -z "$1" ]; then
	echo "Empty sample name supplied to $0, exiting"
	exit 1
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./kraken_weigh_contigs.sh  sample_name  [pre/post] [paired/assembled] run_ID"
	echo "Output is saved to ${processed}/miseq_run_ID_id/sample_name/kraken/(pre/post)assembly/sample_name_kraken_summary_(paired/assembled)"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id supplied to $0, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty kraken version supplied to $0, exiting"
	exit 1
fi

#Sets output folder to the correct path relative to assembly completion
OUTDATADIR="${processed}/${2}/${1}/${3}/postAssembly"
echo "-${OUTDATADIR}-"


#Checks to see if the list file used for calculations exists and exits if it does not
if [[ ! -s "${OUTDATADIR}/${1}_assembled.${3}" ]]; then
	echo "${OUTDATADIR}/${1}_assembled.${3} does not exist"
	exit 1
fi

# Tells script which version of kraken was run (and where to go to get the correct files)
if [[ "${3}" != "kraken" ]] && [[ "${3}" != "kraken2" ]]; then
	echo "Invalid kraken version supplied, exiting"
	exit 4
fi

contig_sizes=()
total_size=0
unclassified=0

sort -t$'\t' -k4,4 -n "${OUTDATADIR}/${1}_assembled.${3}" > "${OUTDATADIR}/${1}_assembled_sorted.${3}"

#Parses the kraken output list line by line
while IFS= read -r line  || [ -n "$line" ]; do
		prefix=$(echo "${line}" | cut -d'	' -f1,2,3,4)
		#echo "${prefix}"
		classified=$(echo "${line}" | cut -d'	' -f1)
		#echo "classified as:${classified}"
		if [[ "${classified}" == "C" ]]; then
			contig_size=$(echo "${line}" | cut -d'	' -f4)
			contig_sizes+=(${contig_size})
			total_size=$(( total_size + contig_size ))
		else
			echo "Contig not classified"
			unclassified=$(( unclassified + 1 ))
		fi
done < "${OUTDATADIR}/${1}_assembled_sorted.${3}"

contig_count=${#contig_sizes[@]}
counter=0
#for contiggy in ${contig_sizes[@]}; do
#	echo "${counter}-${contiggy}"
#	counter=$((counter + 1 ))
#done

smallest=${contig_sizes[0]}
adjusted_contig_count=$(( total_size / smallest ))
echo "Contig count = ${contig_count}"
echo "Total Size = ${total_size}"
echo "unclassified = ${unclassified}"
echo "Smallest contig = ${smallest}"
echo "Adjusted contig count = ${adjusted_contig_count}"

while IFS= read -r line  || [ -n "$line" ]; do
		IFS='	' read -a arr_line <<< "$line"
		#echo "Size:${#arr_line}"
		original_size=${arr_line[3]}
		#echo "3-${arr_line[3]}"
		#echo "OS = ${original_size}, smallest = ${smallest}"
		adjusted_size=$(( original_size / smallest ))
		#echo "Adj_size = ${adjusted_size}"
		arr_line[3]=${adjusted_size}
		echo -e "${arr_line[0]}\t${arr_line[1]}\t${arr_line[2]}\t${arr_line[3]}"

done < "${OUTDATADIR}/${1}_assembled_sorted.${3}"

if [[ ! -s "${OUTDATADIR}/${1}_assembled_weighted.mpa" ]]; then
	echo "${OUTDATADIR}/${1}_assembled_weighted.mpa does not exist, cant do mpa adjustment"
	exit 1
fi

#Script exited gracefully (unless something else inside failed)
exit 0
