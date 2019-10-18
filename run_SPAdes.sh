#!/bin/sh -l

#$ -o run_SPAdes.out
#$ -e run_SPAdes.err
#$ -N run_SPAdes
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Runs SPAdes on sample to align reads into best possible assembly
#
# Usage: ./run_spades.sh sample_name   normal/plasmid    run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/Assembly
#
# Modules required: SPAdes/3.13.0
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
	echo "Empty sample name supplied to run_SPAdes.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_SPAdes.sh sample_name   [normal/plasmid]   run_ID"
	echo "Output by default is sent to ${processed}/miseq_run_ID/sample_name/Assembly"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty run type supplied to run_SPAdes.sh, Should be either normal or  plasmid. Exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty project id supplied to run_SPAdes.sh, exiting"
	exit 1
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${3}/${1}"

##### Non singularity way
###	spades.py --careful --memory "${spades_max_memory}" --only-assembler --pe1-1 "${OUTDATADIR}/trimmed/${1}_R1_001.paired.fq" --pe1-2 "${OUTDATADIR}/trimmed/${1}_R2_001.paired.fq" --pe1-s "${OUTDATADIR}/trimmed/${1}.single.fq" -o "${OUTDATADIR}/Assembly" --phred-offset "${phred}" -t "${procs}"
##### Singularity way
#singularity -s exec -B ${OUTDATADIR}/trimmed:/INPUT -B ${OUTDATADIR}/Assembly:/OUTDIR docker://quay.io/biocontainers/spades:3.13.0--0 spades.py --careful --memory "${spades_max_memory}" --only-assembler --pe1-1 /INPUT/${1}_R1_001.paired.fq --pe1-2 /INPUT/${1}_R2_001.paired.fq --pe1-s /INPUT/${1}.single.fq -o /OUTDIR --phred-offset "${phred}" -t "${procs}"
singularity -s exec -B ${OUTDATADIR}/trimmed:/INPUT -B ${OUTDATADIR}/Assembly:/OUTDIR docker://quay.io/biocontainers/spades:3.13.0--0 spades.py --careful --memory 16 --only-assembler --pe1-1 /INPUT/${1}_R1_001.paired.fq --pe1-2 /INPUT/${1}_R2_001.paired.fq --pe1-s /INPUT/${1}.single.fq -o /OUTDIR --phred-offset "${phred}" -t 6

#Script exited gracefully (unless something else inside failed)
exit 0
