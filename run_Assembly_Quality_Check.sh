#!/bin/sh -l

#$ -o run_Assembly_Quality_Check.out
#$ -e run_Assembly_Quality_Check.err
#$ -N run_Assembly_Quality_Check
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Checks the Assembly quality  using Toms tool and QUAST and comparing the output of both
# 	Important stats are # of contigs, assembly length, n%0 and G/C content
#
# Usage: ./run_Assembly_Quality_Check.sh   sample_name   run_ID [-p]
# 	Optional p flag is to run it only on the plasmid assembly, assuming it is in the default location of the config file and unicycled
#
# Output location: default_config.sh_output_location/run_ID/sample_name/Assembly_Stats(_plasFlow)
#
# Modules required: quast/4.3, python2/2.7.15(loaded by quast/4.3)
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml quast/4.3 Python2/2.7.15

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_Assembly_Quality_Check.sh.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_Assembly_Quality_Check.sh   sample_name run_ID"
	echo "Output is saved to ${processed}/miseq_run_ID/sample_name/Assembly_Stats"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id supplied to run_Assembly_Quality_Check.sh, exiting"
	exit 1
elif [[ "${3}" == "-p" ]]; then
	echo "Only running Assembly quality check on plasFlow assembly"
	do_plasFlow_only="true"
else
	do_plasFlow_only="false"
fi

# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${2}/${1}"


echo "Checking Assembly QC with QUAST"
# Run QUAST
# Save current directory and move to output directory because it doesnt know how to redirect output
owd="$(pwd)"
if [[ "${do_plasFlow_only}" != "true" ]]; then
	# Checks for output folder existence and creates creates if not
	if [ ! -d "$OUTDATADIR/Assembly_Stats" ]; then
		echo "Creating $OUTDATADIR/Assembly_Stats"
		mkdir -p "$OUTDATADIR/Assembly_Stats"
	fi
	cd "${OUTDATADIR}/Assembly_Stats"
	# Call QUAST
	python2 "/apps/x86_64/quast/quast-4.3/quast.py" -o "${OUTDATADIR}/Assembly_Stats" "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta"
	mv "${OUTDATADIR}/Assembly_Stats/report.txt" "${OUTDATADIR}/Assembly_Stats/${1}_report.txt"
	mv "${OUTDATADIR}/Assembly_Stats/report.tsv" "${OUTDATADIR}/Assembly_Stats/${1}_report.tsv"
fi
if [[ -s "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta" ]]; then
	if [ ! -d "$OUTDATADIR/Assembly_Stats_plasFlow" ]; then
		echo "Creating $OUTDATADIR/Assembly_Stats_plasFlow"
 		mkdir -p "$OUTDATADIR/Assembly_Stats_plasFlow"
 	fi
 	python2 "/apps/x86_64/quast/quast-4.3/quast.py" -o "${OUTDATADIR}/Assembly_Stats_plasFlow" "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta"
 	mv "${OUTDATADIR}/Assembly_Stats_plasFlow/report.txt" "${OUTDATADIR}/Assembly_Stats_plasFlow/${1}_report.txt"
 	mv "${OUTDATADIR}/Assembly_Stats_plasFlow/report.tsv" "${OUTDATADIR}/Assembly_Stats_plasFlow/${1}_report.tsv"
 else
	echo "No plasFlow assembly (${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta)"
fi

# Return to original directory
cd "${owd}"

ml -quast/4.3 -Python/2.7.15

#Show that the script seemingly completed successfully
exit 0
