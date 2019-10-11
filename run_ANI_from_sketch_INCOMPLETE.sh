#!/bin/sh -l

#$ -o run_ANI.out
#$ -e run_ANI.err
#$ -N run_ANI
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh

#
# Description: Script to calculate the average nucleotide identity of a sample to a sketch file. *****INCOMPLETE***** at the moment
#
# Usage: ./run_ANI_from_sketch.sh sample_name	run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/ANI/
#
# Modules required: Python3/3.5.2, pyani/0.2.7, mashtree/0.29
#
# V1.0
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml Python3/3.5.2 pyani/0.2.7 mashtree/0.29

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_ANI_from_sketch.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_ANI_from_sketch.sh sample_name run_ID"
	echo "Output is saved to in ${processed}/sample_name/ANI"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty database name supplied to run_ANI.sh. Second argument should be a genus found in ${local_DBs}/ANI/  ...Exiting"
	exit 1
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Started ANI at ${start_time}"

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}"

# Checks to see if an ANI folder already exists and creates it if not
if [ ! -d "$OUTDATADIR/ANI" ]; then
	echo "Creating $OUTDATADIR/ANI"
	mkdir -p "$OUTDATADIR/ANI"
fi

# Checks to see if the local DB ANI folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
if [ ! -d "$OUTDATADIR/ANI/localANIDB" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/ANI/localANIDB"
	mkdir -p "$OUTDATADIR/ANI/localANIDB"
else
	rm -r "$OUTDATADIR/ANI/localANIDB"
	mkdir -p "$OUTDATADIR/ANI/localANIDB"
fi

# Gets persons name to use as email during entrez request to identify best matching sample
me=$(whoami)

# Using the premade sketch file of all refseq assemblies
mash dist "${local_DBs}/aniDB/refseq.genomes.k21s1000.msh" "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" > "${OUTDATADIR}/ANI/${1}_all_refSeq.dists"
sort -k3 -n -o "${OUTDATADIR}/ANI/${1}_all_sorted_refSeq.dists" "${OUTDATADIR}/ANI/${1}_all_refSeq.dists"

exit

counter=0
threshold=0
cp "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" "${OUTDATADIR}/ANI/localANIDB/sample_${1}.fasta"
> "${OUTDATADIR}/ANI/twenty_closest_mash.list"
while IFS='' read -r line || [ -n "$line" ]; do
	source_path=$(echo "${line}" | cut -d'	' -f1)
	source=$(echo "${source_path}" | rev | cut -d'/' -f1 | rev)
	source=${source:0:-4}
	echo ${source}
	distance=$(echo "${line}" | cut -d'	' -f3)
	if [[ ${counter} -gt ${max_ani_samples} ]] && [[ ${distance} != ${threshold} ]]; then
		break
	else
		echo "Going to copy ${source} to localANIDB"
		cp "${source_path}" "${OUTDATADIR}/ANI/localANIDB/${source}.fasta"
		echo "${OUTDATADIR}/ANI/localANIDB/${source}.fasta" >> "${OUTDATADIR}/ANI/twenty_closest_mash.list"
	fi
	counter=$(( counter + 1))
	threshold="${distance}"
done < "${OUTDATADIR}/ANI/${1}_all_sorted.dists"

# Checks for a previous copy of the aniM folder, removes it if found
if [ -d "${OUTDATADIR}/ANI/aniM" ]; then  #checks for and removes old results folder for ANIm
	echo "Removing old ANIm results in ${OUTDATADIR}/ANI/aniM"
	rm -rf "${OUTDATADIR}/ANI/aniM"
fi
exit
#Calls pyani on local db folder
python -V
echo "${OUTDATADIR}/ANI/localANIDB/"
"${shareScript}/pyani_1/average_nucleotide_identity.py" -i "${OUTDATADIR}/ANI/localANIDB/" -o "${OUTDATADIR}/ANI/aniM5" -l "${OUTDATADIR}/ANI/aniM5.log" -m ANIm
#./fastANI --refList "${OUTDATADIR}/ANI/twenty_closest_mash.list" --queryList "${OUTDATADIR}/ANI/twenty_closest_mash.list" -t "${procs}" --matrix -o "${OUTDATADIR}/ANI/matrix_identity.tsv"





This needs to be replaced with the pyani version, or replace the fastANI script before this will work
#./fastANI --refList "${OUTDATADIR}/ANI/twenty_closest_mash.list" --query "${OUTDATADIR}/ANI/localANIDB/sample_${1}.fasta" -t "${procs}"








#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line || [ -n "$line" ]; do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab" ]]; then
	firstline=$(head -n 1 "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab")
else
	echo "No "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab" file, exiting"
	exit 1
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
	definition=$(head -1 "${OUTDATADIR}/ANI/localANIDB/${samples[i]}.fasta")
	# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
	echo "${percents[i+1]} ${samples[i]} ${definition}" >> "${OUTDATADIR}/ANI/best_hits.txt"
done

#Sorts the list in the file based on %id (best to worst)
sort -nr -t' ' -k1 -o "${OUTDATADIR}/ANI/best_hits_ordered.txt" "${OUTDATADIR}/ANI/best_hits.txt"
#Extracts the first line of the file (best hit)
best=$(head -n 1 "${OUTDATADIR}/ANI/best_hits_ordered.txt")
#Creates an array from the best hit
IFS=' ' read -r -a def_array <<< "${best}"
#echo -${def_array[@]}+
#Captures the assembly file name that the best hit came from
best_file=${def_array[1]}
#Formats the %id to standard percentage (xx.xx%)
best_percent=$(awk -v per="${def_array[0]}" 'BEGIN{printf "%.2f", per * 100}')
echo "${best_file}"
# If the best match comes from the additional file, extract the taxonomy from that file
if [[ "${best_file}" = *"_scaffolds_trimmed" ]]; then
	best_outbreak_match=$(echo "${best_file}" | rev | cut -d'_' -f3- | rev)
	while IFS= read -r var  || [ -n "$var" ]; do
		sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
		if [[ "${sample_name}" = "${best_outbreak_match}" ]]; then
			project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
			while IFS= read -r pstats_line || [ -n "$pstats_line" ]; do
					tool=$(echo "${pstats_line}" | cut -d':' -f1 | tr -s " ")
					#echo ":${tool}:"
					if [[ "${tool}" = "weighted Classify " ]]; then
						best_organism_guess=$(echo "${pstats_line}" | cut -d':' -f3 | cut -d' ' -f3,4)
						break 2
					fi
				done < ${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt
		fi
	done < ${5}
# if the best hit comes from the aniDB then pull taxonomy from name
else
	#Extracts the accession number from the definition line
	accession=$(echo "${def_array[2]}" | cut -d' ' -f1  | cut -d'>' -f2)
	#Looks up the NCBI genus species from the accession number
	#best_organism_guess=$(python "${shareScript}/entrez_get_taxon_from_accession.py" -a "${accession}" -e "${me}")

fi

#Creates a line at the top of the file to show the best match in an easily readable format that matches the style on the MMB_Seq log
echo -e "${best_percent}%-${best_organism_guess}(${best_file}.fna)\\n$(cat "${OUTDATADIR}/ANI/best_hits_ordered.txt")" > "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_All).txt"

#Removes the transient hit files
if [ -s "${OUTDATADIR}/ANI/best_hits.txt" ]; then
	rm "${OUTDATADIR}/ANI/best_hits.txt"
#	echo "1"
fi
if [ -s "${OUTDATADIR}/ANI/best_hits_ordered.txt" ]; then
	rm "${OUTDATADIR}/ANI/best_hits_ordered.txt"
#	echo "2"
fi

end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "ENDed ANI at ${end_time}"

#Script exited gracefully (unless something else inside failed)
exit 0
