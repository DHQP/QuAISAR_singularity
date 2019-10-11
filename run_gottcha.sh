#!/bin/sh -l

#$ -o run_gottcha.out
#$ -e run_gottcha.err
#$ -N run_gottcha
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Runs the gottcha classification tool (now only at species level) which identifies the most likely taxonomic classification for the sample
#
# Usage: ./run_gottcha.sh sample_name run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/gottcha/
#
# Modules required: gottcha, perl/5.12.3
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml perl/5.12.3  gottcha

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_gottcha.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_gottcha.sh   sample_name    run_ID"
	echo "Output is saved to ${processed}/miseq_run_ID/sample_name/gottcha/"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id name supplied to run_gottcha.sh, exiting"
	exit 1
fi

# Sets the output folder of gottcha classifier to the gottcha folder under the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}"

# Create necessary output directories
echo "Running gottcha Taxonomic Classifier"
if [ ! -d "$OUTDATADIR/gottcha" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/gottcha"
	mkdir -p "$OUTDATADIR/gottcha"
fi
if [ ! -d "$OUTDATADIR/gottcha/gottcha_S" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/gottcha/gottcha_S"
	mkdir -p "$OUTDATADIR/gottcha/gottcha_S"
fi

if [[ ! -f "${OUTDATADIR}/trimmed/${1}.paired.fq" ]]; then
	if [[ -f "${OUTDATADIR}/trimmed/${1}.paired.fq.gz" ]]; then
		gunzip -c "${OUTDATADIR}/trimmed/${1}.paired.fq.gz" > "${OUTDATADIR}/trimmed/${1}.paired.fq"
	else
		if [[ -f "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" ]]; then
			:
		elif [[ -f "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq.gz" ]]; then
			gunzip -c "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq.gz" > "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq"
		else
			echo "R1_001.paired.fq does not exist (either zipped or not)"
			exit
		fi
		if [[ -f "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" ]]; then
			:
		elif [[ -f "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq.gz" ]]; then
			gunzip -c "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq.gz" > "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq"
		else
			echo "R2_001.paired.fq does not exist (either zipped or not)"
			exit
		fi
		cat "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" > "${OUTDATADIR}/trimmed/${1}.paired.fq"
		if [[ -s "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq.gz" ]]; then
			rm "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq"
		fi
		if [[ -s "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq.gz" ]]; then
			rm "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq.gz"
		fi
	fi
fi


### Gottcha Taxonomy Classifier ### in species mode
gottcha.pl --mode all --outdir "${OUTDATADIR}/gottcha/gottcha_S" --input "${OUTDATADIR}/trimmed/${1}.paired.fq" --database "${gottcha_db}"

gzip "${OUTDATADIR}/trimmed/${1}.paired.fq"

# Return perl to 5.22.1
#. "${shareScript}/module_changers/perl_5123_to_5221.sh"

# Create the krona graphs from each of the analyses
ml krona

ktImportText "${OUTDATADIR}/gottcha/gottcha_S/${1}_temp/${1}.lineage.tsv" -o "${OUTDATADIR}/gottcha/${1}_species.krona.html"

ml -krona
#Create a best hit from gottcha1 file
"${shareScript}/best_hit_from_gottcha1.sh" "${1}" "${2}"

ml -perl/5.12.3 -gottcha

#Script exited gracefully (unless something else inside failed)
exit 0
